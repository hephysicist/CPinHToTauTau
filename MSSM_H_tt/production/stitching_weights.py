from __future__ import annotations
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask
import json

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl     = maybe_import("correctionlib")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

# ---------- small utilities ----------

def _collect_sample_keys_from_cset_json(cset_json_str: str):
  """
  Parse the correctionlib JSON (schemav2) to list all available sample keys.
  Structure: category(era) -> category(sample) -> binning(...)
  Returns: (all_samples: set[str], samples_per_era: dict[era]->set[str])
  """
  try:
    data = json.loads(cset_json_str)
  except Exception:
    return set(), {}
  all_samples = set()
  per_era = {}
  # navigate schemav2
  for corr in data.get("corrections", []):
    if corr.get("name") != None and corr.get("data", {}).get("nodetype") == "category" and corr.get("data", {}).get("input") == "era":
      for era_node in corr["data"].get("content", []):
        era = era_node.get("key")
        per_era.setdefault(era, set())
        inner = era_node.get("value", {})
        if inner.get("nodetype") == "category" and inner.get("input") == "sample":
          for samp_node in inner.get("content", []):
            skey = samp_node.get("key")
            if skey:
              all_samples.add(skey)
              per_era[era].add(skey)
  return all_samples, per_era

def _match_sample_key(dataset_name: str, valid_keys: set[str]) -> str | None:
  """
  Try to map the dataset name to one of the correction's sample keys.
  Strategy: pick the *longest* valid key that is a substring of dataset_name.
  """
  if not dataset_name or not valid_keys:
    return None
  matches = [k for k in valid_keys if k in dataset_name]
  if not matches:
    return None
  # prefer the most specific key (longest)
  return sorted(matches, key=len, reverse=True)[0]

def _choose_var_array(sample_key: str | None, events: ak.Array) -> ak.Array:
  """
  Choose the right LHE variable based on the sample family.
  - W*  -> LHE_Njets
  - DY* -> LHE_NpNLO
  Falls back to whichever is present; if none present, returns zeros with a warning.
  """
  # explicit family check first
  prefer_njets = bool(sample_key and sample_key.startswith("WtoLNu"))
  prefer_npnlo = bool(sample_key and sample_key.startswith("DYto2L"))

  col_njets = "LHE.Njets"
  col_npnlo = "LHE.NpNLO"

  has_njets = has_ak_column(events, col_njets)
  has_npnlo = has_ak_column(events, col_npnlo)

  if prefer_njets and has_njets:
    return ak.values_astype(events.LHE_Njets, np.float32)
  if prefer_npnlo and has_npnlo:
    return ak.values_astype(events.LHE.NpNLO, np.float32)

# ------------------ PRODUCER ------------------

@producer(
  uses={
    "event",
    "LHE.Njets",
    "LHE.NpNLO",       
  },
  produces={
    "stitching_weights",
  },
  mc_only=True,
)
def stitching_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
  """
  Evaluate correctionlib weight:
    weight = corrector(era, sample_key, var)
  where var is LHE_Njets (W) or LHE_NpNLO (DY).
  """
  # figure out era (try config hooks; no hard-coded defaults here)
  era = None
  try:
    era = getattr(self.config_inst.x.stitching_settings, "era", None)
  except Exception:
    pass
  if era is None:
    raise RuntimeError("stitching_weight: 'era' is not set in config (e.g. config.x.stitching_settings).")

  # dataset/sample name: try explicit override first, then derive from dataset name
  sample_key = None
  if sample_key is None:
    ds_name = getattr(self.dataset_inst, "name", None) or getattr(self.dataset_inst, "dataset", None) or ""
    sample_key = _match_sample_key(ds_name, self._stitch_valid_samples)

  if sample_key is None:
    raise RuntimeError(
      "stitching_weight: could not determine sample key. "
      "Either set config.x.stitching_settings.sample explicitly, "
      "or ensure the dataset name contains one of the known samples."
    )

  # choose the right per-event variable and evaluate
  var = _choose_var_array(sample_key, events)
  # correctionlib supports broadcasting: scalars for (era, sample) and array for var
  w = self.stitching.evaluate(era, sample_key, var)

  # enforce float32 and attach
  w = ak.values_astype(w, np.float32)
  events = set_ak_column_f32(events, "stitching_weights", w)
  return events

@stitching_weight.requires
def stitching_weight_requires(self: Producer, reqs: dict) -> None:
  if "external_files" in reqs:
    return
  from columnflow.tasks.external import BundleExternalFiles
  reqs["external_files"] = BundleExternalFiles.req(self.task)

@stitching_weight.setup
def stitching_weight_setup(
  self: Producer,
  reqs: dict,
  inputs: dict,
  reader_targets: InsertableDict,
) -> None:
  bundle = reqs["external_files"]

  import correctionlib
  try:
    # make __call__ behave like .evaluate in older/newer correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
  except Exception:
    pass

  # --- robust load: accept both .json and .json.gz ---
  tgt = bundle.files.stitching_correction
  raw = tgt.load()  # don't pass formatter="gzip" — let law pick by suffix, or return raw bytes/str

  # Normalize to JSON string
  if isinstance(raw, (bytes, bytearray)):
    # detect gzip magic bytes
    if len(raw) >= 2 and raw[0] == 0x1F and raw[1] == 0x8B:
      import gzip
      cjson = gzip.decompress(raw).decode("utf-8")
    else:
      cjson = raw.decode("utf-8")
  elif isinstance(raw, str):
    cjson = raw
  else:
    # in case a formatter returned a Python object
    cjson = json.dumps(raw)

  # cache valid sample keys for name-matching
  self._stitch_valid_samples, self._stitch_samples_per_era = _collect_sample_keys_from_cset_json(cjson)

  # init correction handles
  correction_set = correctionlib.CorrectionSet.from_string(cjson)

  corr_name = getattr(self.config_inst.x.stitching_settings, "corrector", "stitch_weight")
  if corr_name not in correction_set.keys():
    raise RuntimeError(
      f"stitching_weight: corrector '{corr_name}' not found in CorrectionSet. "
      f"Available: {list(correction_set.keys())}"
    )
  self.stitching = correction_set[corr_name]
