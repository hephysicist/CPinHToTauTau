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

from httcp.production.PolarimetricA1 import PolarimetricA1

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

## helper function for the PV theta_GJ rotation
def to_3vec(v):
    return ak.zip({"px": v.px, "py": v.py, "pz": v.pz})

def cross(a, b):
    return ak.zip({
        "px": a.py * b.pz - a.pz * b.py,
        "py": a.pz * b.px - a.px * b.pz,
        "pz": a.px * b.py - a.py * b.px
    })

def unit_3vec(v):
    mag = np.sqrt(v.px**2 + v.py**2 + v.pz**2)
    return ak.zip({
        "px": ak.where(mag > 0, v.px / mag, 0),
        "py": ak.where(mag > 0, v.py / mag, 0),
        "pz": ak.where(mag > 0, v.pz / mag, 0),
    })

def dot_3vec(v1, v2):
    return v1.px * v2.px + v1.py * v2.py + v1.pz * v2.pz

def scale_vec(v, s):
    return ak.zip({
        "px": v.px * s,
        "py": v.py * s,
        "pz": v.pz * s,
    })

def add_vec(v1, v2):
    return ak.zip({
        "px": v1.px + v2.px,
        "py": v1.py + v2.py,
        "pz": v1.pz + v2.pz,
    })

def sub_vec(v1, v2):
    return ak.zip({
        "px": v1.px - v2.px,
        "py": v1.py - v2.py,
        "pz": v1.pz - v2.pz,
    })

def cartesian_to_pt_eta_phi_mass(vector_cartesian):
    """
    Convert an Awkward array with fields (x, y, z, t) representing
    a Lorentz vector in Cartesian coords into a PtEtaPhiMLorentzVectorArray.
    """
    px = vector_cartesian.x
    py = vector_cartesian.y
    pz = vector_cartesian.z
    energy = vector_cartesian.t

    pt = np.sqrt(px**2 + py**2)
    p = np.sqrt(px**2 + py**2 + pz**2)

    epsilon = 1e-12
    eta = 0.5 * np.log((p + pz) / (p - pz + epsilon))
    phi = np.arctan2(py, px)

    mass2 = energy**2 - p**2
    mass2 = ak.where(mass2 > 0, mass2, 0)
    mass = np.sqrt(mass2)

    return ak.zip({
        "pt": pt,
        "eta": eta,
        "phi": phi,
        "mass": mass,
    }, with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)



def rotate_to_gj_max(tau_vis, tau_mtt) -> ak.Array:
    """
    Rotate tau 3-momentum to reduce GJ angle to theta_GJ_max if necessary.

    Parameters:
    - tau_vis: ak.Array of Lorentz vectors (visible decay products)
    - tau_mtt: ak.Array of Lorentz vectors (tau from FastMTT or similar)

    Returns:
    - ak.Array of possibly rotated tau Lorentz vectors (same energy, new direction) if rotation was needed.
    - If no rotation was needed, returns the original tau_mtt.
    """

    # Get the 4-vectors of the tau decay products
    tau_vis_p4 = get_lep_p4(tau_vis)
    tau_mtt_p4 = get_lep_p4(tau_mtt)

    # Unit vectors
    unit = lambda v: ak.where(v.mag > 0, v/v.mag, v/1.)
    tau_vis_dir = unit(tau_vis_p4)
    tau_mtt_dir = unit(tau_mtt_p4)

    # Magnitudes
    tau_vis_mag = tau_vis_p4.mag
    tau_mtt_mag = tau_mtt_p4.mag

    tau_vis_mass = tau_vis_p4.mass
    tau_mtt_mass = tau_mtt_p4.mass

    tau_mtt_energy = tau_mtt_p4.energy

    # Compute cos(theta_GJ)
    cos_theta_gj = tau_vis_dir.dot(tau_mtt_dir)
    cos_theta_gj = np.clip(cos_theta_gj, -1, 1) #ak.where(cos_theta_gj > 1.0, 1.0,ak.where(cos_theta_gj < -1.0, -1.0, cos_theta_gj))
    theta_gj = np.arccos(cos_theta_gj)

    # sin(theta_gj_max)
    sin_theta_gj_max = 0.5 * (tau_mtt_mass**2 - tau_vis_mass**2) / (tau_mtt_mass * tau_vis_mag)
    sin_theta_gj_max = np.clip(sin_theta_gj_max, -1, 1) #ak.where(sin_theta_gj_max > 1.0, 1.0,ak.where(sin_theta_gj_max < -1.0, -1.0, sin_theta_gj_max))
    theta_gj_max = np.arcsin(sin_theta_gj_max)

    # Set default rotated vector
    new_tau_mtt = tau_mtt_p4

    # set coordinate system for rotation
    n1 = to_3vec(tau_vis_dir)
    n3 = unit_3vec(cross(n1, to_3vec(tau_mtt_dir)))
    n2 = unit_3vec(cross(n3, n1))

    P_orth = dot_3vec(to_3vec(tau_mtt_p4), n3)
    P_paral = np.sqrt(tau_mtt_mag**2 - P_orth**2)

    cos_thetaGJmax = np.cos(theta_gj_max)
    sin_thetaGJmax = np.sin(theta_gj_max)

    n_1_paral = unit_3vec(sub_vec(scale_vec(n1, cos_thetaGJmax), scale_vec(n2, sin_thetaGJmax)))
    n_2_paral = unit_3vec(add_vec(scale_vec(n1, cos_thetaGJmax), scale_vec(n2, sin_thetaGJmax)))

    new_dir_1 = unit_3vec(add_vec(scale_vec(n3, P_orth), scale_vec(n_1_paral, P_paral)))
    new_dir_2 = unit_3vec(add_vec(scale_vec(n3, P_orth), scale_vec(n_2_paral, P_paral)))

    # Choose the new direction that aligns better with the original tau_mtt direction
    better_match = dot_3vec(new_dir_1, tau_mtt_dir) > dot_3vec(new_dir_2, tau_mtt_dir)
    new_dir = ak.where(better_match, new_dir_1, new_dir_2)

    new_p3 = scale_vec(new_dir, tau_mtt_mag)

    rotated_tau_mtt = ak.zip({
            'x': ak.from_regular(ak.Array(new_p3.px)),
            'y': ak.from_regular(ak.Array(new_p3.py)),
            'z': ak.from_regular(ak.Array(new_p3.pz)),
            't': ak.from_regular(ak.Array(tau_mtt_energy)),
        }, with_name='PtEtaPhiMLorentzVector', behavior=coffea.nanoevents.methods.vector.behavior)

    # Convert to PtEtaPhiMLorentzVectorArray
    rotated_tau_mtt = cartesian_to_pt_eta_phi_mass(rotated_tau_mtt)

    # Overwrite only where rotates
    mask_theta_gj = theta_gj > theta_gj_max

    rot_tau_mtt = ak.where(mask_theta_gj, rotated_tau_mtt, tau_mtt_p4)

    return rot_tau_mtt


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    tau     = events.hcand_mutau.lep1 
    tau_MTT     = events.hcand_mutau.fastMTT.lep1 #used only in mu_a1_pv
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0    

    p1 = p2 = get_lep_p4(muon)
    r1 = r2 = get_ip_p4(muon)
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
    
    elif pair_decay_ch == "mu_a1_1pr": 
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
       
    elif pair_decay_ch == "mu_a1_3pr":
        charged_pion = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))
        #Create sum of charges of the pions = -1*charge of tau = same charge pion 
        ss_ch, _ = ak.broadcast_arrays(ak.sum(charged_pion.charge, axis=1),
                                         charged_pion.charge)
        ss_pions = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == ss_ch))
        os_pion  = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == -1*ss_ch))
        
        mask_3pr = ((ak.num(ss_pions,axis=1)==2) & (ak.num(os_pion,axis=1)==1))[..., np.newaxis]
        os_pi  = ak.drop_none(ak.mask(get_lep_p4(os_pion),mask_3pr))
        ss_pi1 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:,  :1]),mask_3pr))
        ss_pi2 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, 1:2]),mask_3pr)) 

        m1 = (os_pi + ss_pi1).mass
        m2 = (os_pi + ss_pi2).mass
        
        rho_mass = 0.77526
        dm1 = np.abs(m1 - rho_mass)
        dm2 = np.abs(m2 - rho_mass)
        
        sel_ss_pi = ak.where(dm1 < dm2,
                             ss_pi1,
                             ss_pi2)
        p2 = os_pi #charge of this pion is the same as the charge of tau
        r2 = sel_ss_pi #charge of this pion is opposite to the charge of tau

    elif pair_decay_ch == "mu_a1_3pr_pv":
        charged_pion = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))

        # Create sum of charges of the pions = -1*charge of tau = same charge pion
        ss_ch, _ = ak.broadcast_arrays(ak.sum(charged_pion.charge, axis=1), charged_pion.charge)
        ss_pions = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == ss_ch))
        os_pion  = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == -ss_ch))

        # Select best 2 SS pions and 1 OS pion only for valid combinations
        mask_3pr = ((ak.num(ss_pions,axis=1)==2) & (ak.num(os_pion,axis=1)==1))[..., np.newaxis]
        os_pi  = ak.drop_none(ak.mask(get_lep_p4(os_pion),mask_3pr))
        ss_pi1 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:,  :1]),mask_3pr))
        ss_pi2 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, 1:2]),mask_3pr))

        tau_p4 = rotate_to_gj_max(tau, tau_MTT)

        # Drop events with empty inputs before boosting
        valid_mask = ((ak.num(os_pi) == 1) & (ak.num(tau_p4) == 1) &
                    (ak.num(ss_pi1) == 1) & (ak.num(ss_pi2) == 1))[..., np.newaxis]
        os_pi   = ak.drop_none(ak.mask(os_pi,   valid_mask))
        tau_p4  = ak.drop_none(ak.mask(tau_p4,  valid_mask))
        ss_pi1  = ak.drop_none(ak.mask(ss_pi1,  valid_mask))
        ss_pi2  = ak.drop_none(ak.mask(ss_pi2,  valid_mask))

        vecs_p4 = {
            'p1': os_pi,
            'p2': tau_p4,
            'ss1': ss_pi1,
            'ss2': ss_pi2 # here !  tau_p4 was there before
            }

        zmf_vars = make_boost(vecs_p4)

        a1_pv = PolarimetricA1(zmf_vars['P2'],
                           zmf_vars['P1'],
                           zmf_vars['SS1'],
                           zmf_vars['SS2'],
                           ak.firsts(tau.charge))

        r2 = a1_pv.PVC().boostCM_of_p4(-(vecs_p4['p1'] + vecs_p4['p2']))
        p2 = tau_p4

    final_mask = ((ak.num(p2,axis=1)==1) & (ak.num(r2,axis=1)==1))[...,np.newaxis]   
    p1 = ak.drop_none(ak.mask(p1, final_mask))
    p2 = ak.drop_none(ak.mask(p2, final_mask))
    r1 = ak.drop_none(ak.mask(r1, final_mask))
    r2 = ak.drop_none(ak.mask(r2, final_mask))
    ip2 = ak.drop_none(ak.mask(ip2, final_mask))
    ch1 = ak.drop_none(ak.mask(ch1, final_mask))
    
    if pair_decay_ch == "mu_pi":
        do_phase_shift = ak.zeros_like(r1.energy, dtype=np.bool_) 
    else:
        do_phase_shift = ((p2.energy - r2.energy)/(p2.energy + r2.energy)) < 0
    
    vecs_p4 = {}
    for var in ['p1', 'p2', 'r1', 'r2', 'ip2']:
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



channels = ['mu_pi','mu_rho','mu_a1_1pr','mu_a1_3pr','mu_a1_3pr_pv']

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
