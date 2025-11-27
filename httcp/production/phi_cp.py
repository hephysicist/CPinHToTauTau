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

## helper function for the PV theta_GJ rotation
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


def rotate_to_gj_max(tau_vis, tau_mtt):
    """
    Compute the Gottfried-Jackson angle theta_GJ using 3-vector kinematics.
    Rotate the tau momentum if theta_GJ exceeds the maximal allowed angle theta_max.

    Notes :
    - Angular calculations are performed explicitly without helper functions to correctly handle the array structures.
    - To preserve small-angle precision, vectors are not fully normalised until the final rotation step.
    """

    # --- 4-vectors ---
    vis = get_lep_p4(tau_vis)   # visible decay product
    mtt = get_lep_p4(tau_mtt)   # FastMTT reconstructed tau

    # --- 3-vectors explicitly from pt/eta/phi ---
    vis3 = ak.zip({
        "px": vis.pt * np.cos(vis.phi),
        "py": vis.pt * np.sin(vis.phi),
        "pz": vis.pt * np.sinh(vis.eta),
    })
    mtt3 = ak.zip({
        "px": mtt.pt * np.cos(mtt.phi),
        "py": mtt.pt * np.sin(mtt.phi),
        "pz": mtt.pt * np.sinh(mtt.eta),
    })

    # --- magnitudes ---
    vis_mag = np.sqrt(vis3.px**2 + vis3.py**2 + vis3.pz**2)
    mtt_mag = np.sqrt(mtt3.px**2 + mtt3.py**2 + mtt3.pz**2)

    # --- kinematic theta_GJ via dot product ---
    # Compute theta_GJ preserving small opening angles.
    # Normalizing vectors would erase tiny deviations, giving theta_gj=0 for nearly parallel vectors.
    cos_theta_gj = (vis3.px * mtt3.px + vis3.py * mtt3.py + vis3.pz * mtt3.pz) / (vis_mag * mtt_mag)
    cos_theta_gj = ak.where(cos_theta_gj > 1.0, 1.0,
                            ak.where(cos_theta_gj < -1.0, -1.0, cos_theta_gj))
    theta_gj = np.arccos(cos_theta_gj)

    # --- maximal allowed angle ---
    tau_mass = 1.77686
    vis_m = vis.mass
    denom = 2.0 * tau_mass * vis_mag
    sin_max = ak.where(denom > 0, (tau_mass**2 - vis_m**2) / denom, 0.0)
    sin_max = ak.where(sin_max > 1.0, 1.0, ak.where(sin_max < -1.0, -1.0, sin_max))
    theta_max = np.arcsin(sin_max)

    # --- check which taus need rotation ---
    need_rot = theta_gj > theta_max
    if not ak.any(need_rot):
        return mtt, theta_gj, theta_max, theta_gj

    # --- rotation basis ---
    # Use original 3-vectors to define orthonormal basis
    eps = 1e-9
    n1 = ak.zip({
        "px": ak.where(vis_mag > eps, vis3.px / vis_mag, 0.0),
        "py": ak.where(vis_mag > eps, vis3.py / vis_mag, 0.0),
        "pz": ak.where(vis_mag > eps, vis3.pz / vis_mag, 0.0),
    })
    n3_raw = cross(n1, ak.zip({
        "px": ak.where(mtt_mag > eps, mtt3.px / mtt_mag, 0.0),
        "py": ak.where(mtt_mag > eps, mtt3.py / mtt_mag, 0.0),
        "pz": ak.where(mtt_mag > eps, mtt3.pz / mtt_mag, 0.0),
    }))
    raw_mag = np.sqrt(n3_raw.px**2 + n3_raw.py**2 + n3_raw.pz**2)
    fallback = cross(n1, ak.zip({"px": ak.ones_like(n1.px),
                                 "py": 0*n1.px,
                                 "pz": 0*n1.px}))
    n3 = ak.where(raw_mag > eps, n3_raw, fallback)
    n3 = ak.zip({
        "px": n3.px / np.sqrt(n3.px**2 + n3.py**2 + n3.pz**2),
        "py": n3.py / np.sqrt(n3.px**2 + n3.py**2 + n3.pz**2),
        "pz": n3.pz / np.sqrt(n3.px**2 + n3.py**2 + n3.pz**2),
    })
    n2 = unit_3vec(cross(n3, n1))


    # --- component decomposition along new basis ---
    Porth = dot_3vec(mtt3, n3)
    mtt_p = mtt_mag
    Ppar = np.sqrt(ak.where(mtt_p**2 - Porth**2 > 0, mtt_p**2 - Porth**2, 0.0))

    cosT = np.cos(theta_max)
    sinT = np.sin(theta_max)
    par = unit_3vec(add_vec(scale_vec(n1, cosT), scale_vec(n2, sinT)))

    new_dir = add_vec(scale_vec(par, Ppar), scale_vec(n3, Porth))
    new_dir = unit_3vec(new_dir)

    # --- rotated 4-vector ---
    rot_cart = ak.zip({
        "x": new_dir.px * mtt_p,
        "y": new_dir.py * mtt_p,
        "z": new_dir.pz * mtt_p,
        "t": mtt.energy,
    })
    rotated = cartesian_to_pt_eta_phi_mass(rot_cart)
    theta_rot = ak.where(need_rot, theta_max, theta_gj)
    rotated_tau = ak.where(need_rot, rotated, mtt)

    return rotated_tau, theta_gj, theta_max, theta_rot


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    tau     = events.hcand_mutau.lep1 
    tau_MTT = events.hcand_mutau.fastMTT_BW.lep1
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0    

    p1 = p2 = get_lep_p4(muon)
    r1 = r2 = get_ip_p4(muon)
    ch1 = muon.charge
    ip2 = get_ip_p4(tau)

    # initialise variables for GJ rotation case
    theta_gj = ak.zeros_like(tau.pt) * EMPTY_FLOAT
    theta_max = ak.zeros_like(tau.pt) * EMPTY_FLOAT
    theta_rot = ak.zeros_like(tau.pt) * EMPTY_FLOAT
    
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
       
    elif pair_decay_ch == "mu_a1_3pr_dp_reco":
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

    elif pair_decay_ch.startswith('mu_a1_3pr') and 'gen' not in pair_decay_ch:

        # Select tau four-momentum based on the reconstruction method
        if pair_decay_ch == "mu_a1_3pr_pv_reco":
            tau_p4 = tau
        elif pair_decay_ch == "mu_a1_3pr_pv_mtt":
            tau_p4 = tau_MTT
        elif pair_decay_ch == "mu_a1_3pr_pv_gef":
            # rotate tau to maximise GJ angle
            tau_p4, theta_gj, theta_max, theta_rot = rotate_to_gj_max(tau, tau_MTT)

        # Select charged pions from the tau decay
        charged_pion = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))

        # Identify same-sign (SS) and opposite-sign (OS) pions
        ss_charge, _ = ak.broadcast_arrays(ak.sum(charged_pion.charge, axis=1), charged_pion.charge)
        ss_pions = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == ss_charge))
        os_pion  = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == -ss_charge))

        # Require exactly 2 SS pions and 1 OS pion
        mask_3pr = ((ak.num(ss_pions, axis=1) == 2) & (ak.num(os_pion, axis=1) == 1))[..., np.newaxis]
        os_pi  = ak.drop_none(ak.mask(get_lep_p4(os_pion), mask_3pr))
        tau_p4 = ak.drop_none(ak.mask(get_lep_p4(tau_p4), mask_3pr))
        ss_pi1 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, :1]), mask_3pr))
        ss_pi2 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, 1:2]), mask_3pr))

        # Drop events with missing inputs
        valid_mask = ((ak.num(os_pi) == 1) & (ak.num(tau_p4) == 1) &
                      (ak.num(ss_pi1) == 1) & (ak.num(ss_pi2) == 1))[..., np.newaxis]
        os_pi   = ak.drop_none(ak.mask(os_pi, valid_mask))
        tau_p4  = ak.drop_none(ak.mask(tau_p4, valid_mask))
        ss_pi1  = ak.drop_none(ak.mask(ss_pi1, valid_mask))
        ss_pi2  = ak.drop_none(ak.mask(ss_pi2, valid_mask))

        # Collect vectors for boosting and GJ rotation
        vecs_p4 = {
            'p1': os_pi,
            'p2': tau_p4,
            'ss1': ss_pi1,
            'ss2': ss_pi2
        }

        # Boost vectors to zero-momentum frame of the system
        zmf_vars = make_boost(vecs_p4)

        # Construct polarimetric A1 object
        a1_pv = PolarimetricA1(zmf_vars['P2'],
                               zmf_vars['P1'],
                               zmf_vars['SS1'],
                               zmf_vars['SS2'],
                               ak.firsts(tau.charge))

        # Reconstruct the tau impact-parameter in boosted frame
        r2 = a1_pv.PVC().boostCM_of_p4(-(vecs_p4['p1'] + vecs_p4['p2']))
        p2 = tau_p4  # final corrected tau vector

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
    return vecs_p4, do_phase_shift, ch1, theta_gj, theta_max, theta_rot

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
    

channels = ['mu_pi','mu_rho','mu_a1_1pr',
            'mu_a1_3pr_dp_reco',
            'mu_a1_3pr_pv_reco','mu_a1_3pr_pv_mtt','mu_a1_3pr_pv_gef']

@producer(
    uses={
        "hcand_*","tau_decay_prods_*"},
    produces={
        f"phi_cp_{the_ch}" for the_ch in channels
    }| {
        'theta_gj_mu_a1_3pr_pv_gef',
        'theta_max_mu_a1_3pr_pv_gef',
        'theta_rot_mu_a1_3pr_pv_gef',
    }
)
def phi_cp(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
   
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        vecs_p4, do_phase_shift, ch1, theta_gj, theta_max, theta_rot = prepare_acop_vecs(events, pair_decay_ch=the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4, do_phase_shift, ch1)
        phi_cp = ak.fill_none(ak.firsts(phi_cp,axis=1), EMPTY_FLOAT)
        print(f'Found {ak.sum(phi_cp==EMPTY_FLOAT)}/{len(phi_cp)} phi_cp values that are EMPTY_FLOAT')
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}",phi_cp)

        if the_ch == "mu_a1_3pr_pv_gef":
            theta_gj = ak.fill_none(ak.firsts(theta_gj,axis=1), EMPTY_FLOAT)
            theta_max = ak.fill_none(ak.firsts(theta_max,axis=1), EMPTY_FLOAT)
            theta_rot = ak.fill_none(ak.firsts(theta_rot,axis=1), EMPTY_FLOAT)
            events = set_ak_column_f32(events, f"theta_gj_{the_ch}",theta_gj)
            events = set_ak_column_f32(events, f"theta_max_{the_ch}",theta_max)
            events = set_ak_column_f32(events, f"theta_rot_{the_ch}",theta_rot)

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