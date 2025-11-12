from __future__ import annotations
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from law.util import InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, optional_column as optional

import json
import re



ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl     = maybe_import("correctionlib")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
  uses={
    "event",
    optional("LHE.Njets"),
    optional("LHE.NpNLO"),       
  },
  produces={
    "stitching_weight",
  },
  mc_only=True,
)
def stitching_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    if self.weight_dict is not None:
        var = getattr(events.LHE, self.var)
        weight = np.asarray(np.ones_like(var, dtype=np.float32))
        for nj, the_weight in self.weight_dict.items():
            mask = (var==int(nj))
            weight[mask] = the_weight
    else:
        weight = np.ones_like(events.event, dtype=np.float32)
    events = set_ak_column_f32(events, "stitching_weight", weight)
    return events

# @stitching_weight.requires
# def stitching_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
#   if "external_files" in reqs:
#     return
#   from columnflow.tasks.external import BundleExternalFiles
#   reqs["external_files"] = BundleExternalFiles.req(task)

# @stitching_weight.setup
# def stitching_weight_setup(
#     self: Producer,
#     task: law.Task,
#     reqs: dict,
#     inputs: dict,
#     reader_targets: InsertableDict,
# ) -> None:
#     bundle = reqs["external_files"]
    
#     input_data = bundle.files.stitching.load(formatter='json')
#     self.weight_dict = None
#     dataset = self.dataset_inst.name 
#     if dataset in self.config_inst.x.stitch_samples:
#         year = self.config_inst.x.year
#         tag = self.config_inst.x.long_tag
#         if dataset.startswith("DYto2L"): corr = input_data[f'{year}{tag}']['DY']
#         if dataset.startswith("WtoLNu"): corr = input_data[f'{year}{tag}']['W']
#         self.weight_dict = corr['weights']
#         self.var = corr['condition_variable'].replace('LHE_','')

@stitching_weight.setup
def stitching_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    def get_nj(s: str) -> int:
      m = re.search(r'_(\d+)J_', s)
      return int(m.group(1)) if m else -1
    
    def get_xsec(self, dataset_name):
      d = self.config_inst.get_dataset(dataset_name)
      ecm = self.config_inst.campaign.ecm
      return self.config_inst.get_process(d.processes.names()[0]).get_xsec(ecm)
    self.weight_dict = None
    dataset = self.dataset_inst.name 
    if dataset in self.config_inst.x.stitch_samples:
        start_str = [s for s in ["DYto2L","WtoLNu"] if dataset.startswith(s)]
        self.var = 'NpNLO' if start_str == "DYto2L" else 'Njets'
        self.config_inst.campaign.ecm
        stitch_samples = {get_nj(key): (get_xsec(self, key), 
                                        self.config_inst.get_dataset(key).n_events)
                          for key in self.config_inst.x.stitch_samples 
                          if key.startswith(*start_str)} 
        (xsec_incl, n_evt_incl) = stitch_samples.pop(-1)
        self.weight_dict = {}
        for nj, (xsec,n_evt) in stitch_samples.items():
            num = (xsec/xsec_incl)*n_evt_incl
            self.weight_dict[nj] = num / (num + n_evt)