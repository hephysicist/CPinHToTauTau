import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, DotDict
from law.util import InsertableDict
from columnflow.columnar_util import set_ak_column
from columnflow.types import Any

import law
ak     = maybe_import("awkward")
np     = maybe_import("numpy")
cl     = maybe_import("correctionlib")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={"GenZ.*"},
    produces={"zpt_weight"}, #, "zpt_weight_up", "zpt_weight_down", "zpt_weight_nunc"},
    mc_only=True,
)
def zpt_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    n = len(events)
    one_f32 = np.ones(n, dtype=np.float32)

    sf_nom = ak.Array(one_f32)
    up_mat = ak.Array(np.empty((n, 0), dtype=np.float32))
    down_mat = ak.Array(np.empty((n, 0), dtype=np.float32))
    nunc = 0

    procs = [p.lower() for p in self.dataset_inst.processes.names()]
    if any("dy" in p for p in procs):
        # pt for formulas using log(x); avoid x<=0
        ptll = ak.where(events.GenZ.pt > 0.0, events.GenZ.pt, 1e-6)
        x = ak.to_numpy(ptll)

        # resolve order and base syst
        blob = (self.dataset_inst.name + " " + " ".join(procs)).lower()
        order = getattr(self.config_inst.x, "DY_pTll_order", None)
        if not order:
            order = "NNLO" if "powheg" in blob else ("NLO" if ("amc" in blob or "amcatnlo" in blob) else ("LO" if ("madgraph" in blob or "mg" in blob) else "NLO"))
        base_syst = getattr(self.config_inst.x, "DY_pTll_syst", "nom")

        # nominal
        sf_nom = ak.Array(self.zpt_corr.evaluate(order, x, base_syst).astype(np.float32))

        # how many variations
        nunc = int(self.zpt_nunc_corr.evaluate(order))

        # build up/down matrices: shape (n_events, nunc)
        ups = [self.zpt_corr.evaluate(order, x, f"up{i}").astype(np.float32) for i in range(1, nunc + 1)]
        downs = [self.zpt_corr.evaluate(order, x, f"down{i}").astype(np.float32) for i in range(1, nunc + 1)]

        up_mat = ak.Array(np.stack(ups, axis=1)) if ups else ak.Array(np.empty((n, 0), dtype=np.float32))
        down_mat = ak.Array(np.stack(downs, axis=1)) if downs else ak.Array(np.empty((n, 0), dtype=np.float32))
    events = set_ak_column(events, "zpt_weight", sf_nom, value_type=np.float32)
    # events = set_ak_column(events, "zpt_weight_up", up_mat, value_type=np.float32)
    # events = set_ak_column(events, "zpt_weight_down", down_mat, value_type=np.float32)
    # events = set_ak_column(events, "zpt_weight_nunc", ak.full_like(events.event, nunc, dtype=np.int32), value_type=np.int32)
    return events


@zpt_weight.requires
def zpt_weight_requires(self: Producer, task: law.Task, reqs: dict[str, DotDict[str, Any]], **kwargs) -> None:
    if "external_files" in reqs:
        return
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@zpt_weight.setup
def zpt_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    inputs: dict[str, Any],
    reader_targets: InsertableDict,
) -> None:
    def _norm_path(p):
        if isinstance(p, (tuple, list)):
            return p[0]
        if isinstance(p, dict):
            return p.get("path") or p.get("file") or p.get("filename") or next(iter(p.values()))
        return p

    raw = self.config_inst.x.external_files.zpt_weight
    full_path = _norm_path(raw)
    if not isinstance(full_path, str):
        raise TypeError(f"zpt_weight must resolve to a string path, got: {type(full_path)} = {full_path}")

    cset = cl.CorrectionSet.from_file(full_path)
    self.zpt_corr = cset["DY_pTll_reweighting"]
    self.zpt_nunc_corr = cset["DY_pTll_reweighting_N_uncertainty"]

