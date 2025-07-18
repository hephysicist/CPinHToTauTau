# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

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
    empty_array = ak.zeros_like(array)[..., :0]
    return ak.where(mask, array, empty_array)


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    tau     = events.hcand_mutau.lep1 
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0
    
    p1 = get_lep_p4(muon)
    r1 = get_ip_p4(muon)
    ch1 = muon.charge
    ip2 = get_ip_p4(tau)
    
    if pair_decay_ch == "mu_pi":
        charged_pions = ak.drop_none(ak.mask(tauprod,pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.firsts(tau.charge))
        ss_pions = ak.drop_none(ak.mask(charged_pions, pi_ch == tau_ch))
        sorted_idx = ak.argsort(ss_pions.pt, ascending=False)
        best_pion = ak.drop_none(ak.firsts(ss_pions[sorted_idx],axis=1)[..., np.newaxis])
        p2 = get_lep_p4(best_pion) # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        r2 = get_ip_p4(tau) # Create 4-vectors of tau impact parameters
        r2 = ak.drop_none(ak.mask(r2, r2.rho2 > 0))
        final_mask = ((ak.num(p2,axis=1)==1) & (ak.num(r2,axis=1)==1))[...,np.newaxis]
        p2 = ak.drop_none(ak.mask(p2, final_mask))
        r2 = ak.drop_none(ak.mask(r2, final_mask))
        # For this channel there is no need to do the phase shift, so this arrray is filled with zeros
        do_phase_shift = ak.zeros_like(r1.energy, dtype=np.bool_) 
    elif pair_decay_ch == "mu_rho":
        charged_pions = ak.drop_none(ak.mask(tauprod,pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.firsts(tau.charge))
        ss_pions = ak.drop_none(ak.mask(charged_pions, pi_ch == tau_ch))
        sorted_idx = ak.argsort(ss_pions.pt, ascending=False)
        best_pion = ak.drop_none(ak.firsts(ss_pions[sorted_idx],axis=1)[..., np.newaxis])
        em_particles = ak.drop_none(ak.mask(tauprod,egamma_mask(tauprod)))
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(best_pion)
        # Create 4-vectors of tau impact parameters
        r2 = get_lep_p4(em_particles).sum(1)[...,np.newaxis]
        r2 = ak.drop_none(ak.mask(r2, r2.rho2 > 0)) #This is needed to remove empty events because .sum() returns 4-vector filled with zeros
        final_mask = ((ak.num(p2,axis=1)==1) & (ak.num(r2,axis=1)==1))[...,np.newaxis]
        
        p2 = ak.drop_none(ak.mask(p2, final_mask))
        r2 = ak.drop_none(ak.mask(r2, final_mask))
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0
    
    p1 = ak.drop_none(ak.mask(p1, final_mask))
    r1 = ak.drop_none(ak.mask(r1, final_mask))
    ip2 = ak.drop_none(ak.mask(ip2, final_mask))
    ch1 = ak.drop_none(ak.mask(ch1, final_mask))
    vecs_p4 = {}
    for var in ['p1', 'p2', 'r1', 'r2','ip2']:
        vecs_p4[var] = eval(var)
    return vecs_p4, do_phase_shift, ch1

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
    } 
    # | {
    #     optional(f"phi_cp_{the_ch}_reg1") for the_ch in channels
    # } | {
    #     optional(f"phi_cp_{the_ch}_reg2") for the_ch in channels
    # } | {"phi_cp_incl"} 
    # | {
    #     optional(f"alpha_{the_ch}") for the_ch in channels
    # }
)
def phi_cp(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
   
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        vecs_p4, do_phase_shift, ch1 = prepare_acop_vecs(events, pair_decay_ch=the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4, do_phase_shift, ch1)
        phi_cp = ak.fill_none(ak.firsts(phi_cp,axis=1), EMPTY_FLOAT)
        print(f'Found {ak.sum(phi_cp==EMPTY_FLOAT)}/{len(phi_cp)} phi_cp values that are EMPTY_FLOAT')
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}",phi_cp)
        #alpha = np.ones_like(events.event)*EMPTY_FLOAT
        # alpha_per_ch = produce_alpha(vecs_p4)
        # alpha_per_ch = ak.fill_none(alpha_per_ch, EMPTY_FLOAT)
        
        # reg1_mask = alpha_per_ch >= np.pi/4.
        # reg2_mask = (alpha_per_ch < np.pi/4.) & (alpha_per_ch >= 0.)
        
        # empty_floats = ak.ones_like(phi_cp)*EMPTY_FLOAT
        # for the_reg in ['reg1', 'reg2']:
        #     var = ak.where(eval(f'{the_reg}_mask'), phi_cp, empty_floats)
        #     events = set_ak_column_f32(events, f"phi_cp_{the_ch}_{the_reg}",  var)
      

        # events = set_ak_column_f32(events, f'alpha_{the_ch}', alpha_per_ch)
        # events = set_ak_column_f32(events, "phi_cp_incl", phi_cp)
       
    return events
