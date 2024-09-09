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
from columnflow.production.cms.mc_weight import mc_weight
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
from httcp.production.PolarimetricA1 import PolarimetricA1 as get_pol_vec_for_a1
from httcp.production.dilepton_features import mT

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
warn = maybe_import("warnings")

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


  #helper function to get higgs candidates for single channel. Higgs candidates from different channel are filled with empty hcand obj.
    #Chanel should be "mutau", "etau", "tautau"


#lambda function to get 4-vector from the particle objects
get_lep_p4 = lambda part: ak.zip({ f"{var}" : part[var] for var in ['pt', 'eta', 'phi','mass']},
                                    with_name="PtEtaPhiMLorentzVector", 
                                    behavior=coffea.nanoevents.methods.vector.behavior)
#lambda function to get 4-vector for impact parameter from the particle objects
#Zeroth component of IP vector is set to zero by definition that can be found here: https://www.mdpi.com/2218-1997/8/5/256
get_ip_p4 = lambda part: ak.zip({f'{var}': part[f'IP{var}']for var in ['x', 'y', 'z']} | {'t' : ak.zeros_like(part.IPx)},
                            with_name="LorentzVector", 
                            behavior=coffea.nanoevents.methods.vector.behavior)
   
def get_single_ch_mt(self: Producer, events: ak.Array, tau_decay_channel: str,  **kwargs):
    
    #Produce mask for different channels and tau decay modes
    #Currently the code works for e-tau and mu-tau channels
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    mask = events.channel_id == self.config_inst.get_channel('mutau').id
    mask = mask & events.single_muon_triggered
    presel_tau = ak.firsts(events.hcand[:,1:2])
    presel_tauprod = ak.firsts(events.hcandprod[:,1:2])
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]
    
    mass = (hcand1 + hcand2).mass
    events = set_ak_column(events, "hcand_invm", mass)
    #from IPython import embed; embed()

    #Select the events with taus and muons having ip_significance > 2
 
    mask = mask & ak.fill_none(ak.firsts(events.hcand[:,0:1].ip_sig >= 2), False)
    
    if tau_decay_channel == "mu_pi":
        mask = mask & ak.fill_none(presel_tau.ip_sig >= 2, False)
        mask = mask & ak.fill_none(presel_tau.decayMode == 0, False) # Get only DM0 events
        mask = mask & ak.fill_none(ak.all(presel_tauprod.tauIdx == 0, axis=1), False) #Check that all decay products comming from tau
        mask = mask & ak.fill_none(ak.num(presel_tauprod.pt, axis=1) == 1, False) #Tau has only a single daughter
    elif tau_decay_channel == "mu_rho":
        mask = mask & ak.fill_none(presel_tau.decayMode == 1, False) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.charge != 0 ,axis=1) == 1, False) #Require only one charged product
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.charge == 0 ,axis=1) == 1, False) #Require only one neutral product
    elif tau_decay_channel == "mu_a1_1pr":
        mask = mask & ak.fill_none(presel_tau.decayMode == 1, False)# Get only DM2 events since tau -> a1 + nu -> rho + pi0 + nu 
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.charge != 0 ,axis=1) == 1, False)#Require only one charged product
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.charge == 0 ,axis=1) == 2, False) #Check if all three decay products are comming from the same tau       
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.tauIdx == 0 ,axis=1) == 3, False)
    elif tau_decay_channel == "a1_3pr":
        mask = mask & presel_tau.decayMode == 10 # Get only DM1 events since tau -> pho+nu -> pi^+pi^0+nu
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.charge != 0 ,axis=1) == 3, False) #Require only one charged product
        mask = mask & ak.fill_none(ak.sum(presel_tauprod.tauIdx == 0 ,axis=1) == 3, False) #Check if all three decay products are comming from the same tau
    
    empty_hcand_arr = ak.zeros_like(events.hcand)[...][...,:0]
    empty_hcandprod_arr = ak.zeros_like(events.hcandprod)[...][...,:0]
    sel_hcand = ak.where(mask, events.hcand, empty_hcand_arr)
    sel_hcandprod = ak.where(mask, events.hcandprod, empty_hcandprod_arr)
    #Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
    #lower case variables are defined in laboratory frame
    # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon 
    p1 = get_lep_p4(sel_hcand[:,0:1])
    r1 = get_ip_p4(sel_hcand[:,0:1]) 
    ch1 = sel_hcand[:,0:1].charge
    if tau_decay_channel == "mu_pi":
        p2 = get_lep_p4(ak.flatten(sel_hcandprod,axis=2)) #flatten is needed since there can be multiple decay products for one tau, and DM0 explicitly states that there is only one identified as a decay product. 
        #Create 4-vectors of tau impact parameters
        r2 = get_ip_p4(sel_hcand[:,1:2])
        do_phase_shift = ak.zeros_like(r1.energy, dtype=np.bool_) #For this channel there is no need to do the phase shift, so this arrray is filled with False statements
    elif tau_decay_channel == "mu_rho":
        #Get the charged pion and neutral separately
        charged_pion_mask = sel_hcandprod.charge != 0
        neutral_pion_mask = sel_hcandprod.charge == 0
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(sel_hcandprod,charged_pion_mask)), axis=2))
        r2 = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(sel_hcandprod,neutral_pion_mask)), axis=2))
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0   
    elif tau_decay_channel == "mu_a1_1pr":
        #Get the charged and neutral pions separately
        charged_pion_mask = sel_hcandprod.charge != 0
        neutral_pion_mask = sel_hcandprod.charge == 0
        # for the tau -> a1 decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral system of two neutral particles
        p2 = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(sel_hcandprod,charged_pion_mask)), axis=2))
        neutral_pion = ak.mask(sel_hcandprod,neutral_pion_mask)[:,1:2]
        neutral_pion1 = get_lep_p4(ak.flatten(ak.drop_none(neutral_pion), axis=2)[:,0:1])
        neutral_pion2 = get_lep_p4(ak.flatten(ak.drop_none(neutral_pion), axis=2)[:,1:2])
        r2 = neutral_pion1.add(neutral_pion2)
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0   
            
    # elif tau_decay_channel == "a1_3pr":
    #     #Get flat array for taus and 2d array for tauprods
    #     flat_tau = ak.firsts(sel_hcand[:,1:2])
    #     tauprod = ak.flatten(sel_hcandprod[:,1:2],axis=2)
    #     #broadcast arrays so it's possible to compare charge of tau with the charge of tau decay products
    #     tau, _ = ak.broadcast_arrays(flat_tau, tauprod)
    #     #Get 4-vectors of pions that have same sign with tau and one pion that has opposite sign
    #     ss_pi_mask  = tau.charge == tauprod.charge
    #     os_pi_mask  = tau.charge == -1*tauprod.charge
        
    #     pv_inputs = {}
    #     #There should be two same size pions and one opposite size pion
    #     empty_hcandprod = events.hcandprod[0,0,:0]
    #     empty_hcand     = events.hcand[0,:0]
    #     ss_pi   = ak.drop_none(ak.mask(tauprod,ss_pi_mask), axis=1)
    #     #Take the first pion with the same charge as tau
    #     pv_inputs['ss_pion1_p4']  = get_lep_p4(ak.fill_none(ss_pi[:,0:1], empty_hcandprod, axis=0))
    #     #Take the second pion with the same charge as tau
    #     pv_inputs['ss_pion2_p4']  = get_lep_p4(ak.fill_none(ss_pi[:,1:2], empty_hcandprod, axis=0))
    #     #Take take opposite charge pion
    #     pv_inputs['os_pion_p4']   = get_lep_p4(ak.fill_none(ak.drop_none(ak.mask(tauprod,os_pi_mask), axis=1),empty_hcandprod, axis=0))
    #     pv_inputs['tau_p4']       = get_lep_p4(ak.fill_none(flat_tau,empty_hcand, axis=0))
    #     pv_inputs['tau_charge']   = ak.fill_none(flat_tau.charge,empty_hcand.charge, axis=0)
            
    #     tau_dp_p4 = pv_inputs['ss_pion1_p4'].add(pv_inputs['ss_pion2_p4']).add(pv_inputs['os_pion_p4'])
            
    vecs_p4 = {}
    for var in ['p1','p2', 'r1','r2']:
        exec(f'vecs_p4["{var}"] = {var}')  
    
    return events, vecs_p4, do_phase_shift, ch1


def get_single_ch_tt(self: Producer, events: ak.Array,  tau_decay_channel: str,  **kwargs):
    
    #Produce mask for different channels and tau decay modes
    #Currently the code works for e-tau and mu-tau channels
    events = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    mask = events.channel_id == self.config_inst.get_channel('tautau').id
    hcand_ = ak.with_name(events.hcand, "PtEtaPhiMLorentzVector")
    hcand1 = hcand_[:,0:1]
    hcand2 = hcand_[:,1:2]
    
    mass = (hcand1 + hcand2).mass
    events = set_ak_column(events, "hcand_invm", mass)
    
    leg1_tau =  ak.firsts(events.hcand[:,0:1])
    leg2_tau =  ak.firsts(events.hcand[:,1:2])
    leg1_prod = ak.firsts(events.hcandprod[:,0:1])
    leg2_prod = ak.firsts(events.hcandprod[:,1:2])
     
    if tau_decay_channel == "rho_rho":
        mask = mask & ak.fill_none(leg1_tau.decayMode == 1, False) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & ak.fill_none(ak.sum(leg1_prod.charge != 0 ,axis=1) == 1, False) #Require only one charged product
        mask = mask & ak.fill_none(ak.sum(leg1_prod.charge == 0 ,axis=1) == 1, False) #Require only one neutral product
        
        mask = mask & ak.fill_none(leg2_tau.decayMode == 1, False)
        mask = mask & ak.fill_none(ak.sum(leg2_prod.charge != 0 ,axis=1) == 1, False) #Require only one charged product
        mask = mask & ak.fill_none(ak.sum(leg2_prod.charge == 0 ,axis=1) == 1, False) #Require only one neutral product
    
    empty_hcand_arr = ak.zeros_like(events.hcand)[...][...,:0]
    empty_hcandprod_arr = ak.zeros_like(events.hcandprod)[...][...,:0]
    sel_hcand = ak.where(mask, events.hcand, empty_hcand_arr)
    sel_hcandprod = ak.where(mask, events.hcandprod, empty_hcandprod_arr)
    #Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
    #lower case variables are defined in laboratory frame
    # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon 
    
    ch1 = sel_hcand[:,0:1].charge
    
    leg1_tau =  sel_hcand[:,0:1]
    leg2_tau =  sel_hcand[:,1:2]
    leg1_prod = sel_hcandprod[:,0:1]
    leg2_prod = sel_hcandprod[:,1:2]
    y = ak.ones_like(leg1_tau.pt)
    vecs_p4 = {}
    if tau_decay_channel == "rho_rho":
        for i in range(1,3):
            #Get the charged pion and neutral separately
            
            tau_obj = eval(f'leg{i}_tau')
            tauprod_obj = eval(f'leg{i}_prod')
            charged_pion_mask = tauprod_obj.charge != 0
            neutral_pion_mask = tauprod_obj.charge == 0
            # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
            p = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(tauprod_obj,charged_pion_mask)), axis=2))
            r = get_lep_p4(ak.flatten(ak.drop_none(ak.mask(tauprod_obj,neutral_pion_mask)), axis=2))
            vecs_p4[f"p{i}"] = p
            vecs_p4[f"r{i}"] = r
            y = y * ((p.energy - r.energy)/(p.energy + r.energy))   
    
    return events, vecs_p4, y < 0, ch1

def make_boost(vecs_p4, boostvec = None):
        #Check if  boost vector is provided, othervise calculated it by getting p1 and p2 from the vecs_p4 dict
        if boostvec == None: boostvec_ = vecs_p4['p1'] + vecs_p4['p2']
        else: boostvec_ = boostvec
        #Create a dictionary to store boosted variables (they are defined with upper case names)
        zmf_vars = {}
        for var in vecs_p4.keys():
            exec(f'zmf_vars["{var.upper()}"] = vecs_p4["{var}"].boost(boostvec_.boostvec.negative())')  
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
        
def produce_dy_alpha(
    vecs_p4,
    ch1,
) -> ak.Array:
    empty_p2 = ak.zeros_like(vecs_p4['p2'].pvec.unit)[...][...,:0]
    empty_r2 = ak.zeros_like(vecs_p4['r2'].pvec.unit)[...][...,:0]
    mask = ak.fill_none(ak.firsts(ch1 < 0),False)
    P_ = ak.where(mask, vecs_p4['p1'].pvec.unit, vecs_p4['p2'].pvec.unit)
    R_ = ak.where(mask, vecs_p4['r1'].pvec.unit, vecs_p4['r2'].pvec.unit)
    z = ak.zeros_like(P_)
    z['z'] = 1
    vec1 = (z.cross(P_)).unit
    vec2 = (R_.cross(P_)).unit
    dy_alpha = np.arccos(np.absolute(vec1.dot(vec2)))
    dy_alpha = ak.where(np.isnan(dy_alpha),
                        ak.ones_like(dy_alpha)*EMPTY_FLOAT,
                        dy_alpha)
    return dy_alpha
    

channels = ['mu_pi', 'mu_rho', 'mu_a1_1pr', 'rho_rho']

@producer(
    uses={
        "hcand.*",
        "hcandprod.*",
        "*triggered",
        
        optional("GenTau.*"), optional("GenTauProd.*")},
    produces={
        f"phi_cp_{the_ch}" for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg1") for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg2") for the_ch in channels
    } | {
        f"phi_cp_{the_ch}_2bin" for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg1_2bin") for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg2_2bin") for the_ch in channels
    } |{
        optional(f"alpha_{the_ch}") for the_ch in channels
    } | {"hcand_invm"}
    
)
def get_phi_cp(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    #Now we use the functions defined above for phi_cp calculation
    alpha = np.ones_like(events.event)*EMPTY_FLOAT
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        if 'mu' in the_ch:
            events, vecs_p4, do_phase_shift, ch1 = get_single_ch_mt(self, events, the_ch)
        else:
            events, vecs_p4, do_phase_shift, ch1 = get_single_ch_tt(self, events, the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = ak.firsts(get_acop_angle(zmf_vecs_p4,do_phase_shift))
        alpha_per_ch = produce_dy_alpha(vecs_p4, ch1)
        alpha = ak.where(ak.num(alpha_per_ch)>0,
                                ak.fill_none(ak.firsts(alpha_per_ch),EMPTY_FLOAT),
                                alpha)
        reg1_mask = ak.fill_none(ak.firsts(alpha_per_ch >= np.pi/4.), False)
        reg2_mask = ak.fill_none(ak.firsts((alpha_per_ch <  np.pi/4.) & (alpha_per_ch >= 0.)), False)
        
        phi_cp_2bin = (phi_cp + np.pi/2.) % (2*np.pi)
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}_2bin",  ak.fill_none(phi_cp_2bin, EMPTY_FLOAT))
        
        for the_reg in ['reg1','reg2']:
            var         = ak.where(eval(f'{the_reg}_mask'), phi_cp, EMPTY_FLOAT)
            var_2bin    = ak.where(eval(f'{the_reg}_mask'), phi_cp_2bin, EMPTY_FLOAT)
            
            events = set_ak_column_f32(events, f"phi_cp_{the_ch}_{the_reg}",  var)
            events = set_ak_column_f32(events, f"phi_cp_{the_ch}_{the_reg}_2bin",  var_2bin)
            
        events = set_ak_column_f32(events, f'alpha_{the_ch}', ak.fill_none(alpha, EMPTY_FLOAT))
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}",  ak.fill_none(phi_cp, EMPTY_FLOAT))    
            
    return events
    

@producer(
    uses={
        "hcand.*",
        "hcandprod.*",
        "channel_id",
        optional("GenTau.*"), optional("GenTauProd.*"),
    },
    produces={'channel_id', "hcand",}
)
def mask_spoiled_events(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    # Check if each event contains only one hcand pair
    mask = ak.num(events.hcand.pt, axis=1) == (ak.ones_like(events.event) * 2) 
    if ak.any(~mask):
        warn.warn("hcand array contain multiple candidates for a single channel! \n Masking these events...")
    type_dict = {"pt"            : "float64",
                    "eta"           : "float64",
                    "phi"           : "float64",
                    "mass"          : "float64",
                    "charge"        : "int64",
                    "decayMode"     : "int64",
                    "decayModePNet" : "int64",
                    "rawIdx"        : "int64",
                    "IPx"           : "float64",
                    "IPy"           : "float64",
                    "IPz"           : "float64",
                    "ip_sig"        : "float64"
                    }
    masked_hcand = {}
    masked_channel_id =  ak.values_astype(ak.where(mask, events.channel_id, 255), np.uint8)
    events = set_ak_column(events, 'channel_id', masked_channel_id)
    for var in events.hcand.fields:
        masked_hcand[var] = ak.values_astype(
            ak.fill_none(
                ak.pad_none(
                    events.hcand[var], 2, axis=1),
                1),
            type_dict[var])
        
    hcand_obj = ak.zip(masked_hcand)
    events = set_ak_column(events, 'hcand', hcand_obj)
           
    return events
   
   
    

@producer(
    uses={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        get_phi_cp,
        mask_spoiled_events,
        mc_weight,
        #tauspinner_weight,
        category_ids,
        mT
    },
    produces={
        normalization_weights,
        split_dy,
        pu_weight,
        muon_weight,
        tau_weight,
        get_phi_cp,
        mask_spoiled_events,
        mc_weight,
        #tauspinner_weight,
        category_ids,
        mT
    },
)
def main(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    events = self[mask_spoiled_events](events, **kwargs)
    events = self[mT](events, **kwargs) 
    events = self[category_ids](events, **kwargs) 
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)
        events = self[normalization_weights](events, **kwargs)
        processes = self.dataset_inst.processes.names()
        #if ak.any(['dy' in proc for proc in processes]):
        #    print("Splitting Drell-Yan dataset...")
        #    events = self[split_dy](events, **kwargs)
        events = self[pu_weight](events, **kwargs)
        events = self[muon_weight](events,do_syst=True, **kwargs)
        events = self[tau_weight](events, do_syst=True, **kwargs)
        
        #events = self[tauspinner_weight](events, **kwargs)
    #Create H candidate variables: hcand_invm, hcand_dr
    events = self[get_phi_cp](events, **kwargs)       
    
    return events
