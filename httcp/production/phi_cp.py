# coding: utf-8

"""
Column production methods related to higher-level features.
"""
import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column, optional_column as optional
from httcp.util import get_lep_p4, get_ip_p4
from httcp.production.PolarimetricA1 import PolarimetricA1
from httcp.production.gef import unit, rotate_to_gj_max 

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
warn = maybe_import("warnings")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

def pi0_and_prods_mask(tauprod): return ((np.abs(tauprod.pdgId) == 11) + (tauprod.pdgId == 22) + (tauprod.pdgId == 111) )

def pion_mask(tauprod): return np.abs(tauprod.pdgId) == 211

def get_single_part(array: ak.Array, idx: int) -> ak.Array:
    return ak.firsts(array[:,idx:idx+1])


def prepare_acop_vecs(events: ak.Array, pair_decay_ch):
    PVBS        = events.PVBS
    if pair_decay_ch.endswith("_gen"):
        #from IPython import embed; embed()
        tau = events.gen_lep.lep1
        tauprod = events.gentau_decay_prods_mutau_lep1 
        muon = events.gen_lep.lep0
    else:
        tau = events.hcand_mutau.lep1
        tauprod = events.tau_decay_prods_mutau_lep1
        muon = events.hcand_mutau.lep0
        if (pair_decay_ch == "mu_a1_3pr_pv_mtt") or (pair_decay_ch == "mu_a1_3pr_pv_gef"):
            tauMTT = events.hcand_mutau.fastMTT_BW.lep1
        
    
    ch1 = muon.charge # Take the same charge for reco and gen
    ip2 = get_ip_p4(events.hcand_mutau.lep1) #Take always the IP of reco level tau
    
    r1  = r2 = get_ip_p4(muon)
    p1 = p2 = get_lep_p4(muon)

    # initialise variables for GJ rotation case
    theta_gj, theta_max, theta_rot, empty = (
        ak.zeros_like(tau.pt) * EMPTY_FLOAT for _ in range(4))
    
    dsvpv = ak.zip({"x": empty, "y": empty, "z": empty},
        with_name="Vector3D", behavior=coffea.nanoevents.methods.vector.behavior)
    
    if pair_decay_ch.startswith("mu_pi"):
        charged_pions = ak.drop_none(ak.mask(tauprod,pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.fill_none(ak.firsts(tau.charge),0))
       
        ss_pions = ak.drop_none(ak.mask(charged_pions, pi_ch == tau_ch))
        sorted_idx = ak.argsort(ss_pions.pt, ascending=False)
        best_pion = ss_pions[sorted_idx[:,:1]]
        p2 = get_lep_p4(best_pion) # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        # Create 4-vectors of tau impact parameters
        r2 = get_ip_p4(tau)
        r2 = ak.drop_none(ak.mask(r2, r2.rho2 > 0))
    elif pair_decay_ch.startswith("mu_rho"):
        charged_pions = ak.drop_none(ak.mask(tauprod,pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.fill_none(ak.firsts(tau.charge),0))
        ss_pions = ak.drop_none(ak.mask(charged_pions, pi_ch == tau_ch))
        sorted_idx = ak.argsort(ss_pions.pt, ascending=False)
        best_pion = ss_pions[sorted_idx[:,:1]]
        pi0_and_prods = ak.drop_none(ak.mask(tauprod,pi0_and_prods_mask(tauprod)))
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(best_pion)
        # Create 4-vectors of tau impact parameters
        r2 = get_lep_p4(pi0_and_prods).sum(1)[...,np.newaxis]
        r2 = ak.drop_none(ak.mask(r2, r2.rho2 > 0)) #This is needed to remove empty events because .sum() returns 4-vector filled with zeros
    
    elif pair_decay_ch.startswith("mu_a1_1pr"):
        charged_pions = ak.drop_none(ak.mask(tauprod,pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.fill_none(ak.firsts(tau.charge),0))
        ss_pions = ak.drop_none(ak.mask(charged_pions, pi_ch == tau_ch))
        sorted_idx = ak.argsort(ss_pions.pt, ascending=False)
        best_pion = ss_pions[sorted_idx[:,:1]]
        pi0_and_prods = ak.drop_none(ak.mask(tauprod,pi0_and_prods_mask(tauprod)))
        # for the tau -> rho decay, p1 - is 4-vector of the charged pion and r1 is 4-vector of the neutral pion
        p2 = get_lep_p4(best_pion)
        # Create 4-vectors of tau impact parameters
        r2 = get_lep_p4(pi0_and_prods).sum(1)[...,np.newaxis]
        r2 = ak.drop_none(ak.mask(r2, r2.rho2 > 0)) #This is needed to remove empty events because .sum() returns 4-vector filled with zeros
       
    elif pair_decay_ch.startswith("mu_a1_3pr_dp"):
        charged_pions = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge , ak.fill_none(ak.firsts(tau.charge),0))
        ss_pions = ak.drop_none(ak.mask(charged_pions, charged_pions.charge == tau_ch))
        os_pion  = ak.drop_none(ak.mask(charged_pions, charged_pions.charge == -1*tau_ch))
        
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

    elif pair_decay_ch.startswith('mu_a1_3pr_pv'):
        # Select tau four-momentum based on the reconstruction method
        if (pair_decay_ch == "mu_a1_3pr_pv_mtt") or (pair_decay_ch == "mu_a1_3pr_pv_gef"):
            tau_p4 = get_lep_p4(tauMTT)
        else:    
            tau_p4 = get_lep_p4(tau)
      
        # Identify same-sign (SS) and opposite-sign (OS) pions
        charged_pions = ak.drop_none(ak.mask(tauprod, pion_mask(tauprod)))
        pi_ch, tau_ch = ak.broadcast_arrays(charged_pions.charge ,ak.fill_none(ak.firsts(tau.charge),0))
        ss_pions = ak.drop_none(ak.mask(charged_pions, charged_pions.charge == tau_ch))
        os_pion  = ak.drop_none(ak.mask(charged_pions, charged_pions.charge == -1*tau_ch))

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
        
        tau_charge = ak.drop_none(ak.mask(ak.firsts(tau.charge), mask_3pr))
        
        #Drop events with missing inputs ## same mask as above
        valid_mask = ((ak.num(os_pi) == 1) & (ak.num(tau_p4) == 1) &
                        (ak.num(ss_pi1) == 1) & (ak.num(ss_pi2) == 1))[..., np.newaxis]
        os_pi   = ak.drop_none(ak.mask(os_pi, valid_mask))
        tau_p4  = ak.drop_none(ak.mask(tau_p4, valid_mask))
        ss_pi1  = ak.drop_none(ak.mask(ss_pi1, valid_mask))
        ss_pi2  = ak.drop_none(ak.mask(ss_pi2, valid_mask))
        tau_charge = ak.drop_none(ak.mask(ak.fill_none(ak.firsts(tau.charge),0), valid_mask))

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
            'tau': tau_p4.boost(-beta_tau),
            'os1': os_pi.boost(-beta_tau),
            'ss1': ss_pi1.boost(-beta_tau),
            'ss2': ss_pi2.boost(-beta_tau)}
        a1_pol = PolarimetricA1(vecs_pv['tau'], vecs_pv['os1'], vecs_pv["ss1"], vecs_pv["ss2"], tau_charge)

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

    if pair_decay_ch.startswith("mu_pi"):
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
    

@producer(
    uses={
        "hcand_*","tau_decay_prods_*","PVBS.*",
        optional("gen_lep.*"), optional("gentau_decay_prods_mutau_lep1.*")},
    produces={
        f"phi_cp_*"
    } | {
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
    channels = ['mu_pi',
                'mu_rho',
                'mu_a1_1pr',
                'mu_a1_3pr_dp',
                'mu_a1_3pr_pv']
    if self.dataset_inst.is_mc:
        channels = channels + [the_ch+'_gen' for the_ch in channels]
    channels += ['mu_a1_3pr_pv_mtt','mu_a1_3pr_pv_gef']
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

            empty_event_mask = None
            theta_clean = {}
            theta_vars = {"gj":  theta_gj, "max": theta_max, "rot": theta_rot,}

            for name, theta in theta_vars.items():

                empty_mask = ~np.isfinite(theta)
                print(f"Found {ak.sum(empty_mask)}/{len(theta)} "
                    f"theta_{name} values were non-finite and replaced with EMPTY_FLOAT")

                event_mask = ak.any(empty_mask, axis=1)
                empty_event_mask = (event_mask if empty_event_mask is None
                                    else empty_event_mask | event_mask)
                
                if name == "rot":
                    print("Number of empty events:", ak.sum(empty_event_mask))

                theta_vars[name] = ak.where(np.isfinite(theta), theta, EMPTY_FLOAT)

            theta_gj    = ak.fill_none(ak.firsts(theta_vars["gj"],axis=1), EMPTY_FLOAT)
            theta_max   = ak.fill_none(ak.firsts(theta_vars["max"],axis=1), EMPTY_FLOAT)
            theta_rot   = ak.fill_none(ak.firsts(theta_vars["rot"],axis=1), EMPTY_FLOAT)
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