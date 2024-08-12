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


# @producer(
#     uses={
#         "hcand.*",
#         optional("GenTau.*"), optional("GenTauProd.*"),
#         reArrangeDecayProducts, reArrangeGenDecayProducts,
#         ProduceGenPhiCP, ProducePhiCP,
#     },
#     produces={
#         # new columns
#         "hcand_invm", "hcand_dr",
#         ProduceGenPhiCP, ProducePhiCP,
#     },
# )
# def hcand_features(
#         self: Producer, 
#         events: ak.Array,
#         **kwargs
# ) -> ak.Array:
#     events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
#     hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
#     hcand1 = hcand_[:,0:1]
#     hcand2 = hcand_[:,1:2]
    
#     mass = (hcand1 + hcand2).mass
#     dr = ak.firsts(hcand1.metric_table(hcand2), axis=1)
#     dr = ak.enforce_type(dr, "var * float32")

#     events = set_ak_column(events, "hcand_invm", mass)
#     events = set_ak_column(events, "hcand_dr",   dr)

#     #events, P4_dict     = self[reArrangeDecayProducts](events)
#     #events              = self[ProducePhiCP](events, P4_dict)

#     if "is_signal" in list(self.dataset_inst.aux.keys()):
#         if self.dataset_inst.aux["is_signal"]:
#             pass
#             #events, P4_gen_dict = self[reArrangeGenDecayProducts](events)
#             #events = self[ProduceGenPhiCP](events, P4_gen_dict)
    
#     return events


# def get_single_ch(self: Producer, events: ak.Array, channel: str, tau1_dm: str, tau2_dm=None, **kwargs):
#         mask = events.channel_id == self.config_inst.get_channel(channel).id
#         if tau1_dm == "1pr":
#             mask = mask & ak.flatten(events.hcand[:,1:2].decayMode == 0) # Get only DM0 events
#             mask = mask & ak.flatten(ak.all(events.hcandprod.tauIdx[:,1:2],axis=2) == 0) #Check that all decay products comming from tau
#             mask = mask & ak.flatten(ak.num(events.hcandprod.pt[:,1:2],axis=2) == 1) #Tau has only a single daughter
#         empty_hcand = ak.zeros_like(events.hcand)[...][...,:0]
#         empty_hcandprod = ak.zeros_like(events.hcandprod)[...][...,:0]
#         return ak.where(mask, events.hcand, empty_hcand), ak.where(mask, events.hcandprod, empty_hcandprod)







@producer(
    uses={
        "hcand.*",
        "hcandprod.*",
        optional("GenTau.*"), optional("GenTauProd.*"),
    },
    produces={
        f"phi_cp_mu_{the_ch}" for the_ch in ['pi', 'rho']
    },
)
def get_phi_cp(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    
    #helper function to get higgs candidates for single channel. Higgs candidates from different channel are filled with empty hcand obj.
    #Chanel should be "mutau", "etau", "tautau"
    
    def get_single_ch(self: Producer, events: ak.Array, channel: str, tau_decay_channel: str,  **kwargs):
        
        get_lep_p4 = lambda part: ak.zip({ f"{var}" : part[var] for var in ['pt', 'eta', 'phi','mass']},
                                      with_name="PtEtaPhiMLorentzVector", 
                                      behavior=coffea.nanoevents.methods.vector.behavior)
        get_ip_p4 = lambda part: ak.zip({f'{var}': part[f'IP{var}']for var in ['x', 'y', 'z']} | {'t' : ak.zeros_like(part.IPx)},
                                   with_name="LorentzVector", 
                                   behavior=coffea.nanoevents.methods.vector.behavior)
        
        empty_hcand = ak.zeros_like(events.hcand)[...][...,:0]
        empty_hcandprod = ak.zeros_like(events.hcandprod)[...][...,:0]
    
        mask = events.channel_id == self.config_inst.get_channel(channel).id
       
        if tau_decay_channel == "pi":
            mask = mask & ak.flatten(events.hcand[:,1:2].decayMode == 0) # Get only DM0 events
            mask = mask & ak.flatten(ak.all(events.hcandprod.tauIdx[:,1:2],axis=2) == 0) #Check that all decay products comming from tau
            mask = mask & ak.flatten(ak.num(events.hcandprod.pt[:,1:2],axis=2) == 1) #Tau has only a single daughter
        elif tau_decay_channel == "rho":
            mask = mask & ak.flatten(events.hcand[:,1:2].decayMode == 1) # Get only DM1 events since tau -> pho+nu -> pi^+pi^0+nu
            mask = mask & ak.flatten(ak.sum(events.hcandprod.charge[:,1:2] != 0 ,axis=2) == 1) #Require only one charged product
            mask = mask & ak.flatten(ak.sum(events.hcandprod.charge[:,1:2] == 0 ,axis=2) == 1) #Require only one neutral product    
        elif tau_decay_channel == "a1_3pr":
            mask = mask & ak.flatten(events.hcand[:,1:2].decayMode == 10) # Get only DM1 events since tau -> pho+nu -> pi^+pi^0+nu
            mask = mask & ak.flatten(ak.sum(events.hcandprod.charge[:,1:2] != 0 ,axis=2) == 3) #Require only one charged product
        
        sel_hcand = ak.where(mask, events.hcand, empty_hcand)
        sel_hcandprod = ak.where(mask, events.hcandprod, empty_hcandprod)
        #Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
        #lower case variables are defined in laboratory frame
        if tau_decay_channel == "pi":
            p1 = get_lep_p4(sel_hcand[:,0:1]) 
            p2 = get_lep_p4(ak.flatten(sel_hcandprod,axis=2)) #flatten is needed since there can be multiple decay products for one tau, and DM0 explicitly states that there is only one identified as a decay product. 
            #Create 4-vectors of the impact parameters
            r1 = get_ip_p4(sel_hcand[:,0:1])
            r2 = get_ip_p4(sel_hcand[:,1:2])
            do_phase_shift = ak.zeros_like(r1.x, dtype=np.bool_)
        elif tau_decay_channel == "rho":
            p1 = get_lep_p4(sel_hcand[:,0:1]) 
            r1 = get_ip_p4(sel_hcand[:,0:1])
            #Now we need to get the charged pion and neutral separately
            #from IPython import embed; embed()
            charged_pion_mask = sel_hcandprod.charge != 0
            neutral_pion_mask = sel_hcandprod.charge == 0
            p2 = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(sel_hcandprod,charged_pion_mask)), axis=2))
            r2 = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(sel_hcandprod,neutral_pion_mask)), axis=2))
            do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0   
        elif tau_decay_channel == "a1_3pr":
            p1 = get_lep_p4(sel_hcand[:,0:1]) 
            r1 = get_ip_p4(sel_hcand[:,0:1])
            
            #Get 4-vectors of three pions and sum them up to form a1 system
            pi_plus_mask    = sel_hcandprod.charge > 0
            pi_minus_mask   = sel_hcandprod.charge < 0
            pi_plus  = ak.drop_none(ak.mask(sel_hcandprod,pi_plus_mask))
            pi_minus = ak.drop_none(ak.mask(sel_hcandprod,pi_minus_mask))
            from IPython import embed; embed()
            r2 = get_ip_p4(sel_hcand[:,0:1])
               
        vecs_p4 = {}
        for var in ['p1','p2', 'r1','r2']:
            exec(f'vecs_p4["{var}"] = {var}')  
        
        return vecs_p4, do_phase_shift

    def make_boost(vecs_p4):
        #Create higgs system for ZMF calculation
        h_system = vecs_p4['p1'] + vecs_p4['p2'] 
        #Create a dictionary to store boosted variables (they are defined with upper case names)
        zmf_vars = {}
        for var in ['P1','P2','R1','R2']:
            exec(f'zmf_vars["{var}"] = vecs_p4["{var.lower()}"].boost(h_system.boostvec.negative())')  
        return zmf_vars
    
    def get_acop_angle(vecs_p4, do_phase_shift):
        #Create 4 3-vectors from the vecs_p4 dict
        P1 = vecs_p4['P1'].pvec.unit
        P2 = vecs_p4['P2'].pvec.unit
        R1 = vecs_p4['R1'].pvec.unit
        R2 = vecs_p4['R2'].pvec.unit
        R1_tan = R1.add((P1.multiply(R1.dot(P1))).negative())
        R2_tan = R2.add((P2.multiply(R2.dot(P2))).negative())
        
        O = P2.dot(R1_tan.cross(R2_tan))
        
        raw_phi = np.acos((R1_tan.unit).dot(R2_tan.unit))
        phi_cp = ak.where(O > 0 , raw_phi, 2 * np.pi - raw_phi) 
        phi_cp = ak.where(do_phase_shift, phi_cp + np.pi, phi_cp) 
        #Map  angles into [0,2pi) interval
        phi_cp = ak.where(phi_cp > 2.*np.pi, phi_cp - 2.* np.pi, phi_cp) 
        phi_cp = ak.where(phi_cp < 0, phi_cp + 2.* np.pi, phi_cp)
        return phi_cp
        
    #Now we use the functions defined above for phi_cp calculation
    for the_ch in ['pi', 'rho']:
        print(f"Calculating phi_cp for mu {the_ch}")
        vecs_p4, do_phase_shift = get_single_ch(self, events, "mutau", the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4,do_phase_shift)
        events = set_ak_column_f32(events, f"phi_cp_mu_{the_ch}",  phi_cp)
    return events
    
   
   
    

@producer(
    uses={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        get_phi_cp,
        tauspinner_weight,
    },
    produces={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        get_phi_cp,
        tauspinner_weight,
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    #from IPython import embed; embed()
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
    events = self[get_phi_cp](events, **kwargs)       



    return events
