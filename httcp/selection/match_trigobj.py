# coding: utf-8

"""
Exemplary selection methods.
"""

from typing import Optional

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from IPython import embed
from httcp.util import trigger_object_matching

np = maybe_import("numpy")
ak = maybe_import("awkward")

### Function that stores the trigger.id of the objects (e, mu, and tau) ###

def hlt_path_fired(dictionary):
    assert len(dictionary) > 0, "Empty dictionaly, babushcha alert !!! check triggers !!!"

    max_length = 0
    for key in dictionary.keys():
        temp_length = ak.max(ak.num(dictionary[key], axis=1))
        if temp_length > max_length: max_length = temp_length

    hlt_condition = {}
    for key in dictionary.keys():
        hlt_condition[key] = ak.pad_none(dictionary[key], target=max_length)
        hlt_condition[key] = ak.fill_none(hlt_condition[key],-1)[:,:,None]

    hlt_condition_values = list(hlt_condition.values())
    hlt_condition_values_concat = ak.concatenate(hlt_condition_values, axis=-1)
    HLT_path_fired = ak.max(hlt_condition_values_concat, axis=-1)

    return HLT_path_fired 

@selector(
    uses={
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Muon.pt"    , "Muon.eta"    , "Muon.phi"    , "Muon.mass",
        "Tau.pt",      "Tau.eta",      "Tau.phi",      "Tau.mass",
        "TrigObj.pt",  "TrigObj.eta",  "TrigObj.phi",
    },
    produces={"single_electron_triggered", "cross_electron_triggered", 
              "single_muon_triggered",     "cross_muon_triggered",
              "cross_tau_triggered", 
              "matched_triggerID_e", "matched_triggerID_mu", "matched_triggerID_tau"
          },
    exposed=False
)
def match_trigobj(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        electron_indices: ak.Array,
        muon_indices: ak.Array,
        tau_indices: ak.Array,
        domatch: Optional[bool] = False,
        **kwargs
) -> ak.Array:
    
    # prepare vectors for output vectors
    false_mask = (abs(events.event) < 0)
    single_electron_triggered = false_mask
    cross_electron_triggered  = false_mask
    single_muon_triggered = false_mask
    cross_muon_triggered  = false_mask
    cross_tau_triggered  = false_mask

    electron_indices_dummy = electron_indices[:,:0]
    muon_indices_dummy     = muon_indices[:,:0]
    tau_indices_dummy      = tau_indices[:,:0]
    
    # matched_idx_e = electron_indices[:,:0]
    # matched_idx_mu = muon_indices[:,:0]
    # matched_idx_tau = tau_indices[:,:0]
    # matched_triggerID_e = electron_indices[:,:0]
    # matched_triggerID_mu = muon_indices[:,:0]
    # matched_triggerID_tau = tau_indices[:,:0]
    matched_idx_e   = ak.values_astype(-1 * ak.ones_like(electron_indices), np.int32)
    matched_idx_mu  = ak.values_astype(-1 * ak.ones_like(muon_indices), np.int32)
    matched_idx_tau = ak.values_astype(-1 * ak.ones_like(tau_indices), np.int32)
    matched_triggerID_e   = ak.values_astype(-1 * ak.ones_like(electron_indices), np.int64)
    matched_triggerID_mu  = ak.values_astype(-1 * ak.ones_like(muon_indices), np.int64)
    matched_triggerID_tau = ak.values_astype(-1 * ak.ones_like(tau_indices), np.int64)
    
    hlt_path_fired_e   = {}
    hlt_path_fired_mu  = {}   
    hlt_path_fired_tau = {}      
    if domatch:
        # perform each lepton election step separately per trigger
        for trigger, trigger_fired, leg_masks in trigger_results.x.trigger_data:
            #print(f"trigger: {trigger}")
            #print(f"trigger_fired: {trigger_fired}")
            #print(f"leg_masks:  {leg_masks}")
            
            is_single_el = trigger.has_tag("single_e")
            is_cross_el  = trigger.has_tag("cross_e_tau")
            is_single_mu = trigger.has_tag("single_mu")
            is_cross_mu  = trigger.has_tag("cross_mu_tau")
            is_cross_tau = trigger.has_tag("cross_tau_tau")
            
            if is_single_mu or is_cross_mu:
                mu_matches_leg0 = None
                # muon selection
                muons = events.Muon[muon_indices]
                # start per-muon mask with trigger object matching
                if is_single_mu:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 1
                    assert abs(trigger.legs[0].pdg_id) == 13
                    # match leg 0
                    mu_matches_leg0 = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                    single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                elif is_cross_mu:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 2
                    assert abs(trigger.legs[0].pdg_id) == 13
                    # match leg 0
                    mu_matches_leg0 = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                    cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                hlt_path_fired_mu[trigger.hlt_field]= ak.where(mu_matches_leg0, trigger.id,-1)                    

            if is_single_el or is_cross_el:
                el_matches_leg0 = None
                # electron selection
                electrons = events.Electron[electron_indices]
                # start per-muon mask with trigger object matching
                if is_single_el:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 1
                    assert abs(trigger.legs[0].pdg_id) == 11
                    # match leg 0
                    el_matches_leg0 = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                    single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                elif is_cross_el:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 2
                    assert abs(trigger.legs[0].pdg_id) == 11
                    # match leg 0
                    el_matches_leg0 = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                    cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                    # sel_electron_indices = ak.local_index(electrons[el_matches_leg0])
                hlt_path_fired_e[trigger.hlt_field]= ak.where(el_matches_leg0, trigger.id,-1)

            if is_cross_el or is_cross_mu or is_cross_tau:
                tau_matches = None
                # start per-tau mask with trigger object matching per leg
                taus = events.Tau[tau_indices]
                if is_cross_el or is_cross_mu:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) == 2
                    assert abs(trigger.legs[1].pdg_id) == 15
                    # match leg 1
                    tau_matches_leg1 = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                    tau_matches = tau_matches_leg1
                    if is_cross_el:
                        cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                    elif is_cross_mu:
                        cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                elif is_cross_tau:
                    # catch config errors
                    assert trigger.n_legs == len(leg_masks) >= 2
                    assert abs(trigger.legs[0].pdg_id) == 15
                    assert abs(trigger.legs[1].pdg_id) == 15
                    # match both legs
                    tau_matches_leg0 = trigger_object_matching(taus, events.TrigObj[leg_masks[0]])
                    tau_matches_leg1 = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                    cross_mask = ((tau_matches_leg0 | tau_matches_leg1) 
                                  & ak.any(tau_matches_leg0, axis=1) 
                                  & ak.any(tau_matches_leg1, axis=1))
                    tau_matches = cross_mask
                    cross_tau_triggered = ak.where(trigger_fired & is_cross_tau, True, cross_tau_triggered)
                hlt_path_fired_tau[trigger.hlt_field]= ak.where(tau_matches, trigger.id,-1)

        #from IPython import embed; embed()
    
        triggerID_e = hlt_path_fired(hlt_path_fired_e)
        triggerID_mu = hlt_path_fired(hlt_path_fired_mu)
        triggerID_tau = hlt_path_fired(hlt_path_fired_tau)
    
        mask_triggerID_e = ak.fill_none(triggerID_e > 0, False)
        matched_idx_e = electron_indices[mask_triggerID_e]
        matched_triggerID_e = triggerID_e[mask_triggerID_e]
        
        mask_triggerID_mu = ak.fill_none(triggerID_mu > 0, False)
        matched_idx_mu = muon_indices[mask_triggerID_mu]
        matched_triggerID_mu = triggerID_mu[mask_triggerID_mu]
        
        mask_triggerID_tau = ak.fill_none(triggerID_tau > 0, False)
        matched_idx_tau = tau_indices[mask_triggerID_tau]
        matched_triggerID_tau = triggerID_tau[mask_triggerID_tau]

        electron_indices = matched_idx_e
        muon_indices = matched_idx_mu
        tau_indices = matched_idx_tau

    sel_electron_indices = ak.values_astype(electron_indices, np.int32)
    sel_muon_indices = ak.values_astype(muon_indices, np.int32)
    sel_tau_indices = ak.values_astype(tau_indices, np.int32)


    events = set_ak_column(events, "matched_triggerID_e", matched_triggerID_e)
    events = set_ak_column(events, "matched_triggerID_mu", matched_triggerID_mu)
    events = set_ak_column(events, "matched_triggerID_tau", matched_triggerID_tau)

    events = set_ak_column(events, "single_electron_triggered", single_electron_triggered)
    events = set_ak_column(events, "cross_electron_triggered", cross_electron_triggered)
    events = set_ak_column(events, "single_muon_triggered", single_muon_triggered)
    events = set_ak_column(events, "cross_muon_triggered", cross_muon_triggered)
    events = set_ak_column(events, "cross_tau_triggered", cross_tau_triggered)
                
    
    return events, sel_electron_indices, sel_muon_indices, sel_tau_indices 
