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
    for column in self.weight_columns:
        if not ak.any(['tt' in proc for proc in processes]) and column == 'top_pt_weight':
            print("===")
            print(weight)
            print(Route(column).apply(events),column)
            print("Skipping top_pt_weight for:", processes)
            print(weight)
            print("===")
            continue
        else :
            print("======")
            print(weight)
            weight = weight * Route(column).apply(events)
            print(column, Route(column).apply(events),column)
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
        "tauspinner_weight",
        "zpt_weight",
        "top_pt_weight",
        "trigger_weight_mutau_nom",
    }
    self.uses |= self.weight_columns
    
    # declare shifts that the produced event weight depends on
    shift_sources = {
       "ts",
    }
    self.shifts |= set(get_shifts_from_sources(self.config_inst, *shift_sources))
    
    
    
@weight_producer(
    # both used columns and dependent shifts are defined in init below
    # only run on mc
    mc_only=True,
)
def ff_weight_producer(self: WeightProducer, events: ak.Array, **kwargs) -> ak.Array:
    # build the full event weight
    weight = ak.Array(np.ones(len(events), dtype=np.float32))
    for column in self.weight_columns:
        weight = weight * Route(column).apply(events)
    return events, weight


@ff_weight_producer.init
def ff_weight_producer(self: WeightProducer) -> None:
    # store column names referring to weights to multiply
    self.weight_columns = {
        "normalization_weight",
        "mc_weight",
        "pu_weight",
        "muon_weight_nom",
        "tau_weight_nom",
        "electron_weight_nom",
        "tauspinner_weight",
        "zpt_weight",
        "top_pt_weight",
        "trigger_weight_mutau_nom",
    }
    self.uses |= self.weight_columns
    
    # declare shifts that the produced event weight depends on
    shift_sources = {
       "ts",
    }
    self.shifts |= set(get_shifts_from_sources(self.config_inst, *shift_sources))
