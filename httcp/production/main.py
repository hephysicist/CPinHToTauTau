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
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from columnflow.production.util import attach_coffea_behavior
#from httcp.production.PhiCPNeutralPion import PhiCPNPMethod
from httcp.production.ReArrangeHcandProds import reArrangeDecayProducts, reArrangeGenDecayProducts
from httcp.production.PhiCP_Producer import ProduceDetPhiCP, ProduceGenPhiCP

from httcp.production.weights import muon_weight, tau_weight, get_mc_weight, tauspinner_weight, zpt_weight, electron_weight,fake_factors, trigger_weight_mutau
from httcp.production.sample_split import split_dy
from httcp.production.generatorZ import genZ
from httcp.production.met_recoil import gen_boson, met_recoil
from httcp.production.dilepton_features import hcand_fields,hcand_mt

from httcp.production.phi_cp import phi_cp
from httcp.production.aux_columns import jet_pt_def,jets_taggable, number_b_jet, pion_energy_split
from httcp.production.top_pt_weight import top_pt_weight, gen_parton_top

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

@producer(
    uses={
        attach_coffea_behavior,
        normalization_weights,
        split_dy,
        pu_weight,
        trigger_weight_mutau,
        muon_weight,
        tau_weight,
        electron_weight,
        genZ,
        zpt_weight,
        gen_boson,
        met_recoil,
        get_mc_weight,
        #fake_factors,
        hcand_fields,
        hcand_mt,
        tauspinner_weight,
        phi_cp,
        category_ids,
        gen_parton_top,
        top_pt_weight,
        jet_pt_def,
        jets_taggable,
        number_b_jet,
        "Jet.*",
        pion_energy_split,
        },
    produces={
        attach_coffea_behavior,
        normalization_weights,
        split_dy,
        pu_weight,
        trigger_weight_mutau,
        muon_weight,
        get_mc_weight,
        tau_weight,
        electron_weight,
        genZ,
        zpt_weight,
        gen_boson,
        met_recoil,
        fake_factors,
        hcand_fields,
        hcand_mt,
        tauspinner_weight,
        phi_cp,
        category_ids,
        gen_parton_top,
        top_pt_weight,
        jet_pt_def,
        jets_taggable,
        number_b_jet,
        "Jet.jec_no_jec_diff",
        pion_energy_split,
    },
    # whether weight producers should be added and called
    produce_weights=True,
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)
    print("Producing pion energy split...")
    events = self[pion_energy_split](events, **kwargs)
    print("Producing Jet features...")
    events = set_ak_column_f32(events, "Jet.jec_no_jec_diff", (events.Jet.pt - events.Jet.pt_no_jec))
    print("Producing jet variables for plotting...") 
    events = self[jet_pt_def](events, **kwargs)
    events = self[jets_taggable](events, **kwargs)   
    print("Producing Number of b-jets for categorization")
    events = self[number_b_jet](events, **kwargs)
    print("Producing Hcand features...")
    events = self[hcand_fields](events, **kwargs) 
    if self.dataset_inst.is_mc:
        events = self[gen_boson](events, **kwargs)
    events = self[met_recoil](events,**kwargs)
    events = self[hcand_mt](events, **kwargs) 
    events = self[category_ids](events, **kwargs)

    if self.dataset_inst.is_mc:
        events = self[get_mc_weight](events, **kwargs)
        print("Producing Normalization weights...")
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events,**kwargs)

        events = self[genZ](events, **kwargs)
        print("Z pt reweighting...")
        events = self[zpt_weight](events,**kwargs)
       
        print("Producing PU weights...")          
        events = self[pu_weight](events, **kwargs)
        print("Producing MuTau Trigger SFs weights...")          
        events = self[trigger_weight_mutau](events, **kwargs)      
        print("Producing Muon weights...")
        events = self[muon_weight](events,do_syst = True, **kwargs)
        print("Producing Electron weights...")
        events = self[electron_weight](events,do_syst = True, **kwargs)
        print("Producing Tau weights...")
        events = self[tau_weight](events,do_syst = True, **kwargs)
        print("Producing Tauspinner weights...")
        events = self[tauspinner_weight](events, **kwargs)
        print("Producing GenPartonTop...")
        events = self[gen_parton_top](events, **kwargs)
        top_pt_weight_dummy = ak.where(events.GenPartonTop.pt > 500.0, 500.0, events.GenPartonTop.pt)
        top_pt_weight_dummy = ak.ones_like(top_pt_weight_dummy)
        for variation in ("", "_up", "_down"):
            events = set_ak_column(events, f"top_pt_weight{variation}", top_pt_weight_dummy)
        if (dataset_inst := getattr(self, "dataset_inst", None)) and dataset_inst.has_tag("ttbar"):
            print("Producing Top pT weights...")
            events = self[top_pt_weight](events, **kwargs)
    print("Producing Fake Factor weights...")
    events = self[fake_factors](events, **kwargs)
    print("Producing phi_cp...")
    events = self[phi_cp](events, **kwargs)
    return events
        
@producer(
    uses={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        electron_weight,
        genZ,
        zpt_weight,
        gen_boson,
        met_recoil,
        get_mc_weight,
        hcand_fields,
        hcand_mt,
        tauspinner_weight,
        category_ids,
        gen_parton_top,
        top_pt_weight,
        jet_pt_def,
        pion_energy_split,
        },
    produces={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        get_mc_weight,
        tau_weight,
        electron_weight,
        genZ,
        zpt_weight,
        gen_boson,
        met_recoil,
        hcand_fields,
        hcand_mt,
        tauspinner_weight,
        category_ids,
        gen_parton_top,
        top_pt_weight,
        jet_pt_def,
        pion_energy_split,
    },
)
def ff_method(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing pion energy split...")
    events = self[pion_energy_split](events, **kwargs)
    print("Producing jet variables...") 
    events = self[jet_pt_def](events, **kwargs)
    print("Producing Hcand features...")
    events = self[hcand_fields](events, **kwargs)
    if self.dataset_inst.is_mc:
        events = self[gen_boson](events, **kwargs)
    events = self[met_recoil](events,**kwargs)
    events = self[hcand_mt](events, **kwargs)  
    events = self[category_ids](events, **kwargs)
    if self.dataset_inst.is_mc:
        events = self[get_mc_weight](events, **kwargs)
        print("Producing Normalization weights...")
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events,**kwargs)
        events = self[genZ](events, **kwargs)
        print("Z pt reweighting...")
        events = self[zpt_weight](events,**kwargs)
        print("Producing PU weights...")          
        events = self[pu_weight](events, **kwargs)
        print("Producing Muon weights...")
        events = self[muon_weight](events,do_syst = True, **kwargs)
        print("Producing Electron weights...")
        events = self[electron_weight](events,do_syst = True, **kwargs)
        print("Producing Tau weights...")
        events = self[tau_weight](events,do_syst = True, **kwargs)
        print("Producing Tauspinner weights...")
        events = self[tauspinner_weight](events, **kwargs)
        print("Producing GenPartonTop...")
        events = self[gen_parton_top](events, **kwargs)
        top_pt_weight_dummy = ak.where(events.GenPartonTop.pt > 500.0, 500.0, events.GenPartonTop.pt)
        top_pt_weight_dummy = ak.ones_like(top_pt_weight_dummy)
        for variation in ("", "_up", "_down"):
            events = set_ak_column(events, f"top_pt_weight{variation}", top_pt_weight_dummy)
        if self.dataset_inst.has_tag("ttbar"):
            print("Producing Top pT weights...")
            events = self[top_pt_weight](events, **kwargs)
    return events



