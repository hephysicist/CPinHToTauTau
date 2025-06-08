# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from httcp.util import get_lep_p4, get_ip_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
warn = maybe_import("warnings")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

def egamma_mask(tauprod): return ((np.abs(tauprod.pdgId) == 11) + (tauprod.pdgId == 22))

def pion_mask(tauprod): return np.abs(tauprod.pdgId) == 211

def get_single_part(array: ak.Array, idx: int) -> ak.Array:
    return ak.firsts(array[:,idx:idx+1])

def apply_evt_mask(array: ak.Array, mask: ak.Array) -> ak.Array:
    empty_array = ak.zeros_like(array)[...][..., :0]
    return ak.where(mask, array, empty_array)


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    tau     = events.hcand_mutau.lep1 
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0
    mask = ak.ones_like(tau.pt, dtype=np.bool_)
    mask = mask & (muon.ip_sig >= 1)
    mask = mask & (events.hcand_mutau.mass < 80)
    # Create masks for different channels
    if pair_decay_ch == "mu_pi":
        mask = mask & (tau.ip_sig >= 1)
        mask = mask & (tau.decayModePNet == 0) # Get only DM0 events
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 1) #Make sure that event conteins only one pion

    elif pair_decay_ch == "mu_rho":
        mask = mask & (tau.decayModePNet == 1) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 1) # Require only one charged pion
        mask = mask & (ak.sum(egamma_mask(tauprod), axis=1) >= 1) # Require one or more electron or photon
        
    elif pair_decay_ch == "mu_a1_1pr":
        mask = mask & (tau.decayModePNet == 2) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
        mask = mask & (ak.sum(pion_mask(tauprod), axis=1) == 1) # Require only one charged pion
        mask = mask & (ak.sum(egamma_mask(tauprod), axis=1) >= 1) # Require one or more electron or photon
    
    sel_tau = apply_evt_mask(tau, ak.fill_none(ak.firsts(mask),False))
    sel_tauprod = apply_evt_mask(tauprod, ak.fill_none(ak.firsts(mask),False))
    sel_muon = apply_evt_mask(muon, ak.fill_none(ak.firsts(mask),False))
    
    # Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
    # lower case variables are defined in laboratory frame
    # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon
    p1 = get_lep_p4(sel_muon)
    r1 = get_ip_p4(sel_muon)
    ch1 = sel_muon.charge
    ip2 = get_ip_p4(sel_tau)

    if pair_decay_ch == "mu_pi":
        pion    = ak.mask(sel_tauprod, pion_mask(sel_tauprod))
        p2 = get_lep_p4(pion) # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        r2 = get_ip_p4(sel_tau) # Create 4-vectors of tau impact parameters
        # For this channel there is no need to do the phase shift, so this arrray is filled with zeros
        do_phase_shift = ak.zeros_like(r1.energy, dtype=np.bool_) > 0
    elif pair_decay_ch == "mu_rho":
        charged_pion = ak.mask(sel_tauprod, pion_mask(sel_tauprod))
        em_particles = ak.mask(sel_tauprod, egamma_mask(sel_tauprod))
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(charged_pion)
        # Create 4-vectors of tau impact parameters
        r2 = get_lep_p4(em_particles).sum()[...,np.newaxis]
        r2 = ak.mask(r2, r2.rho2 > 0)
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0
    vecs_p4 = {}
    for var in ['p1', 'p2', 'r1', 'r2','ip2']:
        vecs_p4[var] = eval(f'ak.firsts({var}, axis=1)')
    return vecs_p4, do_phase_shift, ak.firsts(ch1, axis=1)


# def get_single_ch_tt(self: Producer, events: ak.Array,  tau_decay_channel: str,  **kwargs):
    
#     # Produce mask for different channels and tau decay modes
#     events  = ak.Array(events, behavior=coffea.nanoevents.methods.nanoaod.behavior)
#     mask    = events.channel_id == self.config_inst.get_channel('tautau').id
#     tau = []
#     tauprod = []
#     for tau_idx in range(2):
#         tau.append(get_single_part(events.hcand, tau_idx)) 
#         tauprod.append(get_single_part(events.hcandprod, tau_idx))
    
    
#     # Create masks for different channels
    
#     if tau_decay_channel == "pi_pi":
#         for tau_idx in range(2):
#             mask = mask & ak.fill_none(tau[tau_idx].decayMode == 0, False)  # Get only DM0 events
#             mask = mask & ak.fill_none(tau[tau_idx].ip_sig >= 1, False)
#             mask = mask & ak.fill_none(ak.sum(pion_mask(tauprod[tau_idx]), axis=1) == 1, False) #Make sure that event conteins only one pion

#     elif tau_decay_channel == "rho_rho":
#         for tau_idx in range(2):
#             mask = mask & ak.fill_none(tau[tau_idx].decayMode == 1, False) # Get only DM1 events since tau -> rho+nu -> pi^+pi^0+nu
#             mask = mask & ak.fill_none(ak.sum(pion_mask(tauprod[tau_idx]), axis=1) == 1, False)# Require only one charged pion
#             mask = mask & ak.fill_none(ak.sum(egamma_mask(tauprod[tau_idx]), axis=1) >= 1, False) # Require one or more electron or photon
#     hcand     = apply_evt_mask(events.hcand,     mask)
#     hcandprod = apply_evt_mask(events.hcandprod, mask)
#     # Defining a preliminary set of parameters for the function to calculate the acoplanarity angle
#     # lower case variables are defined in laboratory frame
#     # Get p1 and r1 that correspond to kinematic 4-vector and impact parameter vector of the muon
#     # Redefine taus and tauprods now from preselected events
#     tau = []
#     tauprod = []
#     for tau_idx in range(2):
#         tau.append(get_single_part(hcand, tau_idx)) 
#         tauprod.append(get_single_part(hcandprod, tau_idx))
    
#     vecs_p4 = {}
#     ch1 = tau[0].charge
#     if tau_decay_channel == "pi_pi":
#         for tau_idx in range(2):
#             charged_pion = get_single_part(ak.mask(tauprod[tau_idx], pion_mask(tauprod[tau_idx])),0)
#             vecs_p4[f"p{tau_idx+1}"] = get_lep_p4(charged_pion)
#             vecs_p4[f"r{tau_idx+1}"] = get_ip_p4(tau[tau_idx])
#         do_phase_shift = ak.zeros_like(vecs_p4["p1"].energy, dtype=np.bool_)
#     elif tau_decay_channel == "rho_rho":
#         y = []
#         for tau_idx in range(2):
#             charged_pion = get_single_part(ak.mask(tauprod[tau_idx], pion_mask(tauprod[tau_idx])),0)
#             em_particles = ak.mask(tauprod[tau_idx], egamma_mask(tauprod[tau_idx]))
#             # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
#             p = get_lep_p4(charged_pion)
#             r = get_lep_p4(em_particles).sum()
#             r = ak.mask(r, r.rho2 > 0)
#             vecs_p4[f"r{tau_idx+1}"] = p
#             vecs_p4[f"p{tau_idx+1}"] = r
#             y.append((p.energy - r.energy)/(p.energy + r.energy))
#         do_phase_shift = (y[0] * y[1]) < 0

#     return events, vecs_p4, do_phase_shift, ch1


def make_boost(vecs_p4):
    # Create a dictionary to store boosted variables (they are defined with upper case names)
    zmf_vars = {}
    boostvec_ =vecs_p4['p1'].add(vecs_p4['p2'])
    for var in vecs_p4.keys():
        zmf_vars[var.upper()] = vecs_p4[var].boostCM_of_p4(boostvec_)
    return zmf_vars


def get_acop_angle(vecs_p4, do_phase_shift, ch1):
    # Create 4 3-vectors from the vecs_p4 dict
    
    v3 = {}
    for var in vecs_p4.keys():
        v3[var] = vecs_p4[var].to_3D()
    unit = lambda v: ak.where(v.mag > 0, v/v.mag, v/1.)
    for var in v3.keys():
        v3[var] = unit(v3[var])
        
    R1_tan = unit(v3['R1'].add((v3['P1'].multiply(v3['R1'].dot(v3['P1']))).negative()))
    R2_tan = unit(v3['R2'].add((v3['P2'].multiply(v3['R2'].dot(v3['P2']))).negative()))
    
    Pminus = ak.where(ch1 < 0, v3['P1'], v3['P2'])
    Rminus = ak.where(ch1 < 0, R1_tan, R2_tan)
    Rplus  = ak.where(ch1 < 0, R2_tan, R1_tan)
    
    O = Pminus.dot(Rplus.cross(Rminus))
    
    raw_phi = np.acos(Rplus.dot(Rminus))
    
    raw_phi = ak.nan_to_none(raw_phi)
    phi_cp = ak.where(O > 0, raw_phi, 2 * np.pi - raw_phi)
    phi_cp = ak.where(do_phase_shift, phi_cp + np.pi, phi_cp)
    # Map  angles into [0,2pi) interval
    phi_cp = ak.where(phi_cp > 2.*np.pi, phi_cp - 2. * np.pi, phi_cp)
    phi_cp = ak.where(phi_cp < 0, phi_cp + 2. * np.pi, phi_cp)
    return phi_cp


def produce_alpha(vecs_p4) -> ak.Array:
    
    v3 = {}
    for var in vecs_p4.keys():
        v3[var] = vecs_p4[var].to_3D()
    unit = lambda v: ak.where(v.mag > 0, v/v.mag, v/1.)
    for var in v3.keys():
        v3[var] = unit(v3[var])

    P = v3['p2']
    R = v3['ip2']
    z = ak.zeros_like(P)
    z['z'] = 1
    vec1 = unit(z.cross(P))
    vec2 = unit(R.cross(P))
    
    alpha = np.acos(np.absolute(vec1.dot(vec2)))
    
    return alpha


channels = ['mu_pi','mu_rho']


@producer(
    uses={
        "hcand_*","tau_decay_prods_*"},
    produces={
        f"phi_cp_{the_ch}" for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg1") for the_ch in channels
    } | {
        optional(f"phi_cp_{the_ch}_reg2") for the_ch in channels
    } | {"phi_cp_incl"} 
    | {
        optional(f"alpha_{the_ch}") for the_ch in channels
    }

)
def phi_cp(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:

    alpha = np.ones_like(events.event)*EMPTY_FLOAT
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        vecs_p4, do_phase_shift, ch1 = prepare_acop_vecs(events, pair_decay_ch=the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4, do_phase_shift, ch1)
        
        
        phi_cp = ak.fill_none(ak.firsts(phi_cp,axis=1), EMPTY_FLOAT)
       
        alpha_per_ch = produce_alpha(vecs_p4)
        alpha_per_ch = ak.fill_none(alpha_per_ch, EMPTY_FLOAT)
        
        reg1_mask = alpha_per_ch >= np.pi/4.
        reg2_mask = (alpha_per_ch < np.pi/4.) & (alpha_per_ch >= 0.)
        
        empty_floats = ak.ones_like(phi_cp)*EMPTY_FLOAT
        for the_reg in ['reg1', 'reg2']:
            var = ak.where(eval(f'{the_reg}_mask'), phi_cp, empty_floats)
            events = set_ak_column_f32(events, f"phi_cp_{the_ch}_{the_reg}",  var)
      

        events = set_ak_column_f32(events, f'alpha_{the_ch}', alpha_per_ch)
        events = set_ak_column_f32(events, "phi_cp_incl", phi_cp)
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}",phi_cp)
    return events
