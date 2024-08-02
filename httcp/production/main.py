# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.categories import category_ids
from columnflow.production.normalization import normalization_weights
from columnflow.production.cms.pileup import pu_weight
from columnflow.production.cms.seeds import deterministic_seeds
#from columnflow.production.cms.muon import muon_weights
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

#from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProducePhiCP, ProduceGenPhiCP

from httcp.production.weights import muon_weight, tau_weight, tauspinner_weight
from httcp.production.sample_split import split_dy

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@producer(
    uses={
        "hcand.*",
        optional("GenTau.*"), optional("GenTauProd.*"),
        reArrangeDecayProducts, reArrangeGenDecayProducts,
        ProduceGenPhiCP, ProducePhiCP,
    },
    produces={
        # new columns
        "hcand_invm", "hcand_dr",
        ProduceGenPhiCP, ProducePhiCP,
    },
)
def hcand_features(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]
    
    mass = (hcand1 + hcand2).mass
    dr = ak.firsts(hcand1.metric_table(hcand2), axis=1)
    dr = ak.enforce_type(dr, "var * float32")

    events = set_ak_column(events, "hcand_invm", mass)
    events = set_ak_column(events, "hcand_dr",   dr)

    events, P4_dict     = self[reArrangeDecayProducts](events)
    events              = self[ProducePhiCP](events, P4_dict)

    if "is_signal" in list(self.dataset_inst.aux.keys()):
        if self.dataset_inst.aux["is_signal"]:
            events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
            events = self[ProduceGenPhiCP](events, P4_gen_dict)
    
    return events


@producer(
    uses={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        hcand_features,
        tauspinner_weight,
    },
    produces={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        hcand_features,
        tauspinner_weight,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    if self.dataset_inst.is_mc:
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        #if ak.any(['dy' in proc for proc in processes]):
        #print("Splitting Drell-Yan dataset...")
        #events = self[split_dy](events, **kwargs)
        events = self[pu_weight](events, **kwargs)
        events = self[muon_weight](events,do_syst=True, **kwargs)
        events = self[tau_weight](events, do_syst=True, **kwargs)
        if 'is_signal' in list(self.dataset_inst.aux.keys()): 
            events = self[tauspinner_weight](events, cp_hypo=self.config_inst.x.cp_hypo, **kwargs)
    #Create H candidate variables: hcand_invm, hcand_dr
    #events = self[hcand_features](events, **kwargs)       



    return events
