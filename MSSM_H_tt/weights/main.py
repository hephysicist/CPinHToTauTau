# coding: utf-8

"""
Example event weight producer.
"""

from columnflow.weight import WeightProducer, weight_producer
from columnflow.util import maybe_import
from columnflow.config_util import get_shifts_from_sources
from columnflow.columnar_util import Route

ak = maybe_import("awkward")
np = maybe_import("numpy")


@weight_producer(
    # both used columns and dependent shifts are defined in init below
    # only run on mc
    mc_only=True,
)
def main(self: WeightProducer, events: ak.Array, **kwargs) -> ak.Array:
    # build the full event weight
    processes = self.dataset_inst.processes.names()
    
    weight = ak.Array(np.ones(len(events), dtype=np.float32))

    _stitch_allow = set(self.config_inst.x.stitch_DYto2L_samples)

    _dataset_name = getattr(self.dataset_inst, "name", "")

    for column in self.weight_columns:
        # keep your existing skip rule for top_pt_weight
        if ((self.dataset_inst.has_tag("ttbar") ^ True) & (column == 'top_pt_weight')):
            print("===")
            print(weight)
            print(Route(column).apply(events), column)
            print("Skipping top_pt_weight for:", _dataset_name)
            print(weight)
            print("===")
            continue

        # NEW: only apply stitching_weights to the selected samples
        if column == 'stitching_weights' and _dataset_name not in _stitch_allow:
            print("===")
            print(weight)
            print(Route(column).apply(events), column)
            print("Skipping stitching_weights for:", _dataset_name)
            print(weight)
            print("===")
            continue
        
        # default: apply the weight
        print("======")
        print(weight)
        weight = weight * Route(column).apply(events)
        print(column, Route(column).apply(events), column)
        print(weight)
        print("======")


    process_id = events.process_id
    Z_ee_weight = 1
    if ak.any(['dy' in proc for proc in processes]) and self.dataset_inst.campaign.x.tag == 'postEE':
        print("Applying an ad hoc weight to Z->ee...")
        weight = ak.where(process_id==51001,weight*Z_ee_weight,weight)
        print("The ad hoc weight to Z->ee was applied succefully...")

    return events, weight



@main.init
def main_init(self: WeightProducer) -> None:
    # store column names referring to weights to multiply
    self.weight_columns = {
        "normalization_weight",
        "mc_weight",
        "pu_weight",
        "muon_weight_nom",
        "tau_weight_nom",
        "electron_weight_nom",
        "zpt_weight",
        "btag_weight_SF_nom",
        "top_pt_weight",
        "Trigger_SF_nom",
        "stitching_weights",
    }
    self.uses |= self.weight_columns
    
    # declare shifts that the produced event weight depends on
    shift_sources = {}
    self.shifts |= set(get_shifts_from_sources(self.config_inst, *shift_sources))