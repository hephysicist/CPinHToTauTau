# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from httcp.util import get_lep_p4, get_ip_p4, to_pt_eta_phi_m
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

# helper functions for the unit vectors as the method is buggy
def unit(v, eps=1e-9):
    v3 = v.to_3D()
    mag = v3.mag

    return ak.zip({
            "x": ak.where(mag > eps, v3.x / mag, 1.0),
            "y": ak.where(mag > eps, v3.y / mag, 1.0),
            "z": ak.where(mag > eps, v3.z / mag, 1.0),},
        with_name="Vector3D", behavior=coffea.nanoevents.methods.vector.behavior)


def rotate_to_gj_max(vis, mtt, dsvpv):
    """
    Enforce the physical Gottfried–Jackson limit for the tau decay system. 
    The function computes theta_GJ between the visible system (vis) 
    and the tau direction approximated by the FastMTT output (mtt).
    If theta_GJ exceeds theta_max (the kinematic upper bound),
    the tau momentum is rotated inside the decay plane such that
    theta_GJ = theta_max.
    
    Notes :
    - All angular and cross-product operations are performed in 3D to maintain numerical stability. 
      Only the tau energy is carried through from the original four-vector for the final reconstruction.
    - To preserve small-angle precision, vectors are not fully normalised until the final rotation step.
    """

    vis3 = vis.to_3D()

    tau_mass = 1.77686
    tau_comi = ak.zip({
            "pt": mtt.pt,
            "eta": dsvpv.eta,
            "phi": dsvpv.phi,
            "mass": tau_mass,},
        with_name="PtEtaPhiMLorentzVector",behavior=coffea.nanoevents.methods.vector.behavior,)
    tau_comi3 = tau_comi.to_3D()

    # --- Compute theta_GJ preserving small opening angles ---
    # Normalising vectors would erase tiny deviations, giving theta_gj=0 for nearly parallel vectors.
    cos_theta_gj = vis3.dot(tau_comi3) / (vis3.mag * tau_comi3.mag)
    cos_theta_gj = ak.where(cos_theta_gj > 1.0, 1.0,
                            ak.where(cos_theta_gj < -1.0, -1.0, cos_theta_gj))
    theta_gj = np.arccos(cos_theta_gj)

    # --- maximal allowed angle : sin(theta_max) = (m_tau^2 - m_vis^2) / (2 m_tau |p_vis|) ---
    denom = 2.0 * tau_mass * vis3.mag
    sin_max = ak.where(denom > 0, (tau_mass**2 - vis.mass**2) / denom, 0.0)
    sin_max = ak.where(sin_max > 1.0, 1.0, ak.where(sin_max < -1.0, -1.0, sin_max))
    theta_max = np.arcsin(sin_max)

    # --- rotation basis ---
    eps = 1e-9
    n1 = unit(vis3)
    tau_comi_unit = unit(tau_comi3)

    n3 = unit(n1.cross(tau_comi_unit))
    n2 = unit(n3.cross(n1))

    # --- decompose mtt momentum in this basis (signed components) ---
    Porth = tau_comi3.dot(n3)
    P1 = tau_comi3.dot(n1)
    P2 = tau_comi3.dot(n2)
    p_plane_sq = P1**2 + P2**2
    Ppar = np.sqrt(ak.where(p_plane_sq > 0, p_plane_sq, 0.0))

    # --- compute cos/sin of theta_max ---
    cosT = np.cos(theta_max)
    sinT = np.sin(theta_max)

    # --- two candidate in-plane directions at theta_max (both signs) ---
    n_par_1 = unit(n1 * cosT - n2 * sinT)
    n_par_2 = unit(n1 * cosT + n2 * sinT)

    # --- reconstruct two candidate full directions keeping Porth unchanged ---
    new_dir_1 = unit(ak.where(Ppar > 0, n3 * Porth + n_par_1 * Ppar, n1))
    new_dir_2 = unit(ak.where(Ppar > 0, n3 * Porth + n_par_2 * Ppar, n1))

    # --- choose the candidate closer to the original tau flight direction ---
    choose_1 = new_dir_1.dot(tau_comi_unit) > new_dir_2.dot(tau_comi_unit)
    new_dir = ak.where(choose_1, new_dir_1, new_dir_2)

    # --- final spatial momentum (preserve original |p|) & rebuild 4-vector ---
    new_p = new_dir * tau_comi3.mag
    E_new = np.sqrt(tau_comi3.mag**2 + tau_mass**2) # equal results to mtt.energy

    tau_rot = ak.zip({
            "px": new_p.x,
            "py": new_p.y,
            "pz": new_p.z,
            "E": E_new},
        with_name="PtEtaPhiMLorentzVector", behavior=coffea.nanoevents.methods.vector.behavior,)
    tau_rot = to_pt_eta_phi_m(tau_rot)
    
    # --- check which taus need rotation ---
    need_rot = theta_gj > theta_max

    tau_final = ak.where(need_rot, tau_rot, tau_comi)
    tau3_final = tau_final.to_3D()

    cos_theta_rot = vis3.dot(tau3_final) / (vis3.mag * tau3_final.mag)
    cos_theta_rot = ak.where(cos_theta_rot > 1.0, 1.0,
                    ak.where(cos_theta_rot < -1.0, -1.0, cos_theta_rot))
    theta_rot = np.arccos(cos_theta_rot)

    return tau_final, theta_gj, theta_max, theta_rot


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    tau     = events.hcand_mutau.lep1
    tau_MTT = events.hcand_mutau.fastMTT_BW.lep1
    tauprod = events.tau_decay_prods_mutau_lep1
    muon    = events.hcand_mutau.lep0   
    PVBS      = events.PVBS

    p1 = p2 = get_lep_p4(muon)
    r1 = r2 = get_ip_p4(muon)
    ch1 = muon.charge
    ip2 = get_ip_p4(tau)

    # initialise variables for GJ rotation case
    theta_gj, theta_max, theta_rot, empty = (
        ak.zeros_like(tau.pt) * EMPTY_FLOAT for _ in range(4))
    
    dsvpv = ak.zip({"x": empty, "y": empty, "z": empty},
        with_name="Vector3D", behavior=coffea.nanoevents.methods.vector.behavior)
    
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
        
        mask_3pr = ((ak.num(ss_pions,axis=1)>=2) & (ak.num(os_pion,axis=1)>=1))[..., np.newaxis]
        ss_sorted_idx = ak.argsort(ss_pions.pt, axis=1, ascending=False)
        ss_pions      = ss_pions[ss_sorted_idx]
        os_sorted_idx = ak.argsort(os_pion.pt, axis=1, ascending=False)
        os_pion       = os_pion[os_sorted_idx]
        os_pi  = ak.drop_none(ak.mask(get_lep_p4(os_pion[:, 0:1]),mask_3pr))
        ss_pi1 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:,  0:1]),mask_3pr))
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

    elif pair_decay_ch.startswith('mu_a1_3pr_pv') and 'gen' not in pair_decay_ch:

        # Select tau four-momentum based on the reconstruction method
        if pair_decay_ch == "mu_a1_3pr_pv_reco":
            tau_p4 = get_lep_p4(tau)
        elif pair_decay_ch == "mu_a1_3pr_pv_mtt":
            tau_p4 = get_lep_p4(tau_MTT)
        elif pair_decay_ch == "mu_a1_3pr_pv_gef":
            tau_p4 = get_lep_p4(tau_MTT)

        # Select charged pions from the tau decay
        charged_pion = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))

        # Identify same-sign (SS) and opposite-sign (OS) pions
        ss_charge, _ = ak.broadcast_arrays(ak.sum(charged_pion.charge, axis=1), charged_pion.charge)
        ss_pions = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == ss_charge))
        os_pion  = ak.drop_none(ak.mask(charged_pion, charged_pion.charge == -ss_charge))

        # Require exactly 2 SS pions and 1 OS pion
        mask_3pr = ((ak.num(ss_pions, axis=1) >= 2) & (ak.num(os_pion, axis=1) >= 1))[..., np.newaxis]
        ss_sorted_idx = ak.argsort(ss_pions.pt, axis=1, ascending=False)
        ss_pions      = ss_pions[ss_sorted_idx]
        os_sorted_idx = ak.argsort(os_pion.pt, axis=1, ascending=False)
        os_pion       = os_pion[os_sorted_idx]
        os_pi  = ak.drop_none(ak.mask(get_lep_p4(os_pion[:, 0:1]), mask_3pr))
        tau_p4 = ak.drop_none(ak.mask(tau_p4, mask_3pr))
        ss_pi1 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, 0:1]), mask_3pr))
        ss_pi2 = ak.drop_none(ak.mask(get_lep_p4(ss_pions[:, 1:2]), mask_3pr))

        # m1 = (os_pi + ss_pi1).mass
        # m2 = (os_pi + ss_pi2).mass
        
        # rho_mass = 0.77526
        # dm1 = np.abs(m1 - rho_mass)
        # dm2 = np.abs(m2 - rho_mass)
        
        # sel_ss_pi = ak.where(dm1 < dm2, ss_pi1, ss_pi2)
        # single_pi = ak.where(dm1 > dm2, ss_pi1, ss_pi2)
        
        # Drop events with missing inputs ## same mask as above
        valid_mask = ((ak.num(os_pi) == 1) & (ak.num(tau_p4) == 1) &
                        (ak.num(ss_pi1) == 1) & (ak.num(ss_pi2) == 1))[..., np.newaxis]
                        #(ak.num(sel_ss_pi) == 1) & (ak.num(single_pi) == 1))[..., np.newaxis]
        os_pi   = ak.drop_none(ak.mask(os_pi, valid_mask))
        tau_p4  = ak.drop_none(ak.mask(tau_p4, valid_mask))
        ss_pi1  = ak.drop_none(ak.mask(ss_pi1, valid_mask))
        ss_pi2  = ak.drop_none(ak.mask(ss_pi2, valid_mask))
        # ss_pi1  = ak.drop_none(ak.mask(sel_ss_pi, valid_mask))
        # ss_pi2  = ak.drop_none(ak.mask(single_pi, valid_mask))
        tau_charge = ak.drop_none(ak.mask(ak.firsts(tau.charge), valid_mask))

        # GJ rotation and boost in Higgs rest frame
        if pair_decay_ch == "mu_a1_3pr_pv_gef":
            tau_vis = os_pi + ss_pi1 + ss_pi2

            # # Boost to Higgs rest frame
            # muon = get_lep_p4(muon)
            # muon = ak.drop_none(ak.mask(muon, valid_mask))
            # higgs_p4 = muon + tau_p4
            # beta_H = higgs_p4.boostvec
            # tau_vis_H = tau_vis.boost(-beta_H)
            # tau_p4_H = tau_p4.boost(-beta_H)

            # # rotate tau to max GJ angle if kinematic limit is exceeded
            # tau_p4_H_rot, theta_gj, theta_max, theta_rot = rotate_to_gj_max(tau_vis_H, tau_p4_H)

            # # boost back to lab frame
            # tau_p4 = tau_p4_H_rot.boost(beta_H)

            tau = ak.drop_none(ak.mask(tau, valid_mask))
            PVBS = ak.drop_none(ak.mask(PVBS, valid_mask))
            dsvpv = ak.zip({
                "x" : tau.refitSVx - PVBS.x,
                "y" : tau.refitSVy - PVBS.y,
                "z" : tau.refitSVz - PVBS.z}, with_name="Vector3D", behavior=coffea.nanoevents.methods.vector.behavior)
            dsvpv = ak.drop_none(ak.mask(dsvpv, valid_mask))
    
            tau_p4, theta_gj, theta_max, theta_rot = rotate_to_gj_max(tau_vis, tau_p4, dsvpv)
        
        # Collect vectors for boosting into tau rest frame
        beta_tau = tau_p4.boostvec
        vecs_pv = {
            'p1': os_pi.boost(-beta_tau),
            'p2': tau_p4.boost(-beta_tau),
            'ss1': ss_pi1.boost(-beta_tau),
            'ss2': ss_pi2.boost(-beta_tau)}

        a1_pol = PolarimetricA1(vecs_pv['p2'], vecs_pv['p1'], vecs_pv["ss1"], vecs_pv["ss2"], tau_charge)

        # Reconstruct the tau impact-parameter in boosted frame
        pv = -a1_pol.PVC().pvec

        pv = ak.zip({
            "x": pv.x, 
            "y": pv.y, 
            "z": pv.z, 
            "t": ak.zeros_like(pv.x)}, with_name="LorentzVector", behavior=coffea.nanoevents.methods.vector.behavior)
        r2 = pv.boost(beta_tau)

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
    return vecs_p4, do_phase_shift, ch1, theta_gj, theta_max, theta_rot, dsvpv

def make_boost(vecs_p4, boostvec_=None):
    # Create a dictionary to store boosted variables (they are defined with upper case names)
    zmf_vars = {}
    if boostvec_ is None:
        boostvec_ = vecs_p4['p1'].add(vecs_p4['p2'])
        for var in vecs_p4.keys():
            zmf_vars[var.upper()] = vecs_p4[var].boostCM_of_p4(boostvec_)
    return zmf_vars


def get_acop_angle(vecs_p4, do_phase_shift, ch1):

    v3 = {k: vecs_p4[k].to_3D() for k in vecs_p4}
    v3 = {k: unit(v) for k, v in v3.items()}

    R1_tan = unit(v3['R1'] - v3['P1'] * v3['R1'].dot(v3['P1']))
    R2_tan = unit(v3['R2'] - v3['P2'] * v3['R2'].dot(v3['P2']))
    
    Pminus = ak.where(ch1 < 0, v3['P1'], v3['P2'])
    Rminus = ak.where(ch1 < 0, R1_tan, R2_tan)
    Rplus  = ak.where(ch1 < 0, R2_tan, R1_tan)
    
    # Compute orientation via scalar triple product: O = Pminus · (Rplus × Rminus)
    # determines the sign of the angle
    O = Pminus.dot(Rplus.cross(Rminus))
    
    # Compute angle between Rplus and Rminus
    # Clip cosine to [-1, 1] to prevent arccos(x) with |x| > 1 due to floating point errors
    cos_angle = Rplus.dot(Rminus)
    cos_angle = ak.where(cos_angle > 1.0, 1.0, ak.where(cos_angle < -1.0, -1.0, cos_angle))
    raw_phi = np.arccos(cos_angle)
    
    # Apply orientation correction: use raw_phi if O > 0, otherwise use 2π - raw_phi
    phi_cp = ak.where(O > 0, raw_phi, 2.0 * np.pi - raw_phi)
    
    # Apply optional phase shift of π
    phi_cp = ak.where(do_phase_shift, phi_cp + np.pi, phi_cp)
    
    # Normalize angle to [0, 2π) interval
    phi_cp = ak.where(phi_cp >= 2.0 * np.pi, phi_cp - 2.0 * np.pi, phi_cp)
    phi_cp = ak.where(phi_cp < 0.0, phi_cp + 2.0 * np.pi, phi_cp)
    
    return phi_cp


def produce_alpha(vecs_p4) -> ak.Array:
    
    v3 = {key: unit(vecs_p4[key]) for key in vecs_p4}

    P = v3['p2']
    R = v3['ip2']
    #z = ak.zeros_like(P)
    z = ak.zip({"x": 0*P.x, "y": 0*P.x, "z": 1 + 0*P.x},
                with_name="Vector3D", behavior=coffea.nanoevents.methods.vector.behavior)
    vec1 = unit(z.cross(P))
    vec2 = unit(R.cross(P))
    
    alpha = np.acos(np.absolute(vec1.dot(vec2)))
    
    return alpha
    

channels = ['mu_pi','mu_rho','mu_a1_1pr',
            'mu_a1_3pr_dp_reco',
            'mu_a1_3pr_pv_reco','mu_a1_3pr_pv_mtt','mu_a1_3pr_pv_gef']

@producer(
    uses={
        "hcand_*","tau_decay_prods_*","PVBS.*"},
    produces={
        f"phi_cp_{the_ch}" for the_ch in channels
    }| {
        'theta_gj_mu_a1_3pr_pv_gef',
        'theta_max_mu_a1_3pr_pv_gef',
        'theta_rot_mu_a1_3pr_pv_gef',
        'dsvpv_x','dsvpv_y','dsvpv_z','dsvpv_mag',
    }
)
def phi_cp(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
   
    for the_ch in channels:
        print(f"Calculating phi_cp for {the_ch}")
        vecs_p4, do_phase_shift, ch1, theta_gj, theta_max, theta_rot, dsvpv = prepare_acop_vecs(events, pair_decay_ch=the_ch)
        zmf_vecs_p4 = make_boost(vecs_p4)
        phi_cp = get_acop_angle(zmf_vecs_p4, do_phase_shift, ch1)
        phi_cp = ak.fill_none(ak.firsts(phi_cp,axis=1), EMPTY_FLOAT)
        print(f'Found {ak.sum(phi_cp==EMPTY_FLOAT)}/{len(phi_cp)} phi_cp values that are EMPTY_FLOAT')
        print(f"Found {ak.sum(~np.isfinite(phi_cp))}/{len(phi_cp)} phi_cp values were non-finite and replaced with EMPTY_FLOAT")
        phi_cp = ak.where(np.isfinite(phi_cp), phi_cp, EMPTY_FLOAT)
        events = set_ak_column_f32(events, f"phi_cp_{the_ch}",phi_cp)

        if the_ch == "mu_a1_3pr_pv_gef":
            theta_gj  = ak.fill_none(ak.firsts(theta_gj,axis=1), EMPTY_FLOAT)
            theta_max = ak.fill_none(ak.firsts(theta_max,axis=1), EMPTY_FLOAT)
            theta_rot = ak.fill_none(ak.firsts(theta_rot,axis=1), EMPTY_FLOAT)
            dsvpv_x     = ak.fill_none(ak.firsts(dsvpv.x,axis=1), EMPTY_FLOAT)
            dsvpv_y     = ak.fill_none(ak.firsts(dsvpv.y,axis=1), EMPTY_FLOAT)
            dsvpv_z     = ak.fill_none(ak.firsts(dsvpv.z,axis=1), EMPTY_FLOAT)
            dsvpv_mag   = ak.fill_none(ak.firsts(dsvpv.mag,axis=1), EMPTY_FLOAT)
            events = set_ak_column_f32(events, f"theta_gj_{the_ch}",theta_gj)
            events = set_ak_column_f32(events, f"theta_max_{the_ch}",theta_max)
            events = set_ak_column_f32(events, f"theta_rot_{the_ch}",theta_rot)
            events = set_ak_column_f32(events, f"dsvpv_x",dsvpv_x)
            events = set_ak_column_f32(events, f"dsvpv_y",dsvpv_y)
            events = set_ak_column_f32(events, f"dsvpv_z",dsvpv_z)
            events = set_ak_column_f32(events, f"dsvpv_mag",dsvpv_mag)

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
