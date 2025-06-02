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

from MSSM_H_tt.production.weights import muon_weight, tau_weight, get_mc_weight, zpt_weight, electron_weight, trigger_sf
from MSSM_H_tt.production.sample_split import split_dy
from MSSM_H_tt.production.generatorZ import generatorZ
from MSSM_H_tt.production.dilepton_features import hcand_fields,hcand_mt

from MSSM_H_tt.production.aux_columns import jet_pt_def,jets_taggable,number_b_jet
from MSSM_H_tt.production.btag_SF import btag_weight_SF
from MSSM_H_tt.production.top_pt_weight import top_pt_weight, gen_parton_top
from MSSM_H_tt.production.D_zeta import D_zeta
from MSSM_H_tt.production.met_recoil_correction import gen_boson, met_recoil
#from MSSM_H_tt.production.DY_recoil_unc import DY_pTll_recoil_unc
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
        muon_weight,
        tau_weight,
        electron_weight,
        trigger_sf,
        generatorZ,
        zpt_weight,
        get_mc_weight,
        hcand_fields,
        hcand_mt,
        category_ids,
        number_b_jet,
        jet_pt_def,
        jets_taggable,
        btag_weight_SF,
        gen_parton_top,
        top_pt_weight,
        D_zeta,
        gen_boson,
        met_recoil,
        trigger_sf,
        #DY_pTll_recoil_unc,
        },
    produces={
        attach_coffea_behavior,
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        get_mc_weight,
        tau_weight,
        electron_weight,
        trigger_sf,
        generatorZ,
        zpt_weight,
        hcand_fields,
        hcand_mt,
        category_ids,
        number_b_jet,
        jet_pt_def,
        jets_taggable,
        btag_weight_SF,
        gen_parton_top,
        top_pt_weight,
        D_zeta,
        gen_boson,
        met_recoil,
        trigger_sf,
        #DY_pTll_recoil_unc,     
    },
    # whether weight producers should be added and called
    produce_weights=True,
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    processes = self.dataset_inst.processes.names()
    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)
    print("Producing jet variables for plotting...") 
    events = self[jet_pt_def](events, **kwargs)
    events = self[jets_taggable](events, **kwargs)   
    print("Producing Number of b-jets for categorization")
    events = self[number_b_jet](events, **kwargs)
    print("Producing D_zeta features...")
    events = self[D_zeta](events, **kwargs)
    print("Producing Hcand features...")
    events = self[hcand_fields](events, **kwargs) 
    events = self[category_ids](events, **kwargs)
    
    if (self.dataset_inst.is_mc & (self.config_inst.channels.names()[0] != 'emu')):
        if ak.any(['dy' in proc for proc in processes]):
            print("Splitting Drell-Yan dataset...")
            events = self[split_dy](events,**kwargs)
    
    
    print("Producing D_zeta features...")
    events = self[D_zeta](events, **kwargs)

    if self.dataset_inst.is_mc:
        #if ak.any(['dy' in proc for proc in processes]) or ak.any(['wj' in proc for proc in processes]):
        print("Applying recoil corrections on DY and W+jets samples...")
        events = self[gen_boson](events, **kwargs)
        events = self[met_recoil](events,**kwargs)   
            # events = self[DY_pTll_recoil](events,**kwargs)
            # print("Evaluate recoil corrections uncertainties for DY and W+jets samples...")
            # events = self[DY_pTll_recoil_unc](events,**kwargs)
        print("Getting mc weights...")
        events = self[get_mc_weight](events, **kwargs)
        print("Producing Normalization weights...")
        events = self[normalization_weights](events, **kwargs)
        events = self[generatorZ](events, **kwargs) 
        print("Z pt reweighting...")
        events = self[zpt_weight](events,**kwargs)
        print("Producing PU weights...")          
        events = self[pu_weight](events, **kwargs)
        print("Producing Muon weights...")
        events = self[muon_weight](events,do_syst = True, **kwargs)
        print("Producing Electron weights...")
        events = self[electron_weight](events,do_syst = True, **kwargs)
        print("Producing SFs from efficiencies...")
        events = self[trigger_sf](events, **kwargs)
        print("Producing Tau weights...")
        events = self[tau_weight](events,do_syst = True, **kwargs)
        print("Producing btag weights...")
        events = self[btag_weight_SF](events,do_syst = True,**kwargs)
        print("Producing GenPartonTop...")
        events = self[gen_parton_top](events, **kwargs)
        top_pt_weight_dummy = ak.where(events.GenPartonTop.pt > 500.0, 500.0, events.GenPartonTop.pt)
        top_pt_weight_dummy = ak.ones_like(top_pt_weight_dummy)
        for variation in ("", "_up", "_down"):
            events = set_ak_column(events, f"top_pt_weight{variation}", top_pt_weight_dummy)
        if (dataset_inst := getattr(self, "dataset_inst", None)) and dataset_inst.has_tag("ttbar"):
            print("Producing Top pT weights...")
            events = self[top_pt_weight](events, **kwargs)
    print("Producing mT distributions...") 
    events = self[hcand_mt](events, **kwargs) 

    return events
    