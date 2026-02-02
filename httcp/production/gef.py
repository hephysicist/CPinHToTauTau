# coding: utf-8

"""
A producer to apply the Global Event Fit (GEF) aka Gottfried–Jackson 
angle rotation to the hadronic lepton candidates in the hcand columns.
This steps proceeds after FastMTT. FastMTT corrects the tau momenta 
amplitude (pT, E), while GEF corrects their direction (eta, phi).
"""

import functools
from columnflow.util import maybe_import
from httcp.util import to_pt_eta_phi_m

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

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
