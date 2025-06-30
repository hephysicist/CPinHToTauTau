# coding: utf-8

"""
Collection of helpers
"""

from __future__ import annotations


import law
import order as od
from typing import Any
from columnflow.util import maybe_import
from columnflow.columnar_util import ArrayFunction, deferred_column
from columnflow.util import DotDict
from columnflow.columnar_util import set_ak_column
from columnflow.production import Producer, producer


np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



from columnflow.types import Any
from columnflow.columnar_util import ArrayFunction, deferred_column


@deferred_column
def IF_NANO_V9(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V11(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 11 else None

# lambda function to get 4-vector of a lepton
def get_lep_p4(part): return ak.zip({f"{var}": part[var] for var in ['pt', 'eta', 'phi', 'mass']},
                                    with_name="PtEtaPhiMLorentzVector",
                                    behavior=coffea.nanoevents.methods.vector.behavior)

def get_p2(part): return ak.zip({f"{var}": part[var] for var in ['pt', 'phi']},
                                    with_name="PolarTwoVector",
                                    behavior=coffea.nanoevents.methods.vector.behavior)

# lambda function to get 4-vector of impact parameter from the particle objects
# Zeroth component of IP vector is set to zero by definition that can be found here: https://www.mdpi.com/2218-1997/8/5/256
def get_ip_p4(part): return ak.zip({f'{var}': part[f'IP{var}']for var in ['x', 'y', 'z']} | {'t': ak.zeros_like(part.IPx)},
                                   with_name="LorentzVector",
                                   behavior=coffea.nanoevents.methods.vector.behavior)# # lambda function to get 4-vector from the particle objects


@deferred_column
def IF_NANO_V12(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 12 else None


@deferred_column
def IF_RUN_2(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.run == 2 else None


@deferred_column
def IF_RUN_3(self: ArrayFunction.DeferredColumn, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.run == 3 else None


@deferred_column
def IF_DATASET_HAS_LHE_WEIGHTS(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return None if func.dataset_inst.has_tag("no_lhe_weights") else self.get()


@deferred_column
def IF_DATASET_IS_DY(
    self: ArrayFunction.DeferredColumn,
    func: ArrayFunction,
) -> Any | set[Any]:
    if getattr(func, "dataset_inst", None) is None:
        return self.get()

    return self.get() if func.dataset_inst.has_tag("is_dy") else None

def transverse_mass(lepton: ak.Array, met: ak.Array) -> ak.Array:
    dphi_lep_met = lepton.delta_phi(met)
    mt = np.sqrt(2 * lepton.pt * met.pt * (1 - np.cos(dphi_lep_met)))
    return mt


def trigger_object_matching(
    vectors1: ak.Array,
    vectors2: ak.Array,
    threshold: float = 0.5,
    axis: int = 2,
) -> ak.Array:
    """
    Helper to check per object in *vectors1* if there is at least one object in *vectors2* that
    leads to a delta R metric below *threshold*. The final reduction is applied over *axis* of the
    resulting metric table containing the full combinatorics. When *return_all_matches* is *True*,
    the matrix with all matching decisions is returned as well.
    """
    # delta_r for all combinations
    dr = vectors1.metric_table(vectors2)
    # check per element in vectors1 if there is at least one matching element in vectors2
    any_match = (dr < threshold)
    return any_match


def hlt_path_fired(events,dictionary):
    HLT_path_fired = -1*ak.ones_like(ak.local_index(events.event), dtype=np.int64)
    if len(dictionary) > 0:
        # Start by collecting all the arrays from the dictionary
        # Collect arrays from the dictionary
        array_list = []
        for key in dictionary.keys():
            array_list.append(dictionary[key])

        # Stack all the arrays along the second axis (axis=1)
        hlt_condition_values_concat = np.stack(array_list, axis=1)
        HLT_path_fired = ak.fill_none(ak.max(hlt_condition_values_concat, axis=-1),-1)
    return HLT_path_fired 
    
# Define the fill_missing function
def fill_missing(arr, fill_value):
    return ak.fill_none(arr, fill_value)
   
def get_dataset_lfns(
        dataset_inst: od.Dataset,
        shift_inst: od.Shift,
        dataset_key: str,
) -> list[str]:
    # destructure dataset_key into parts and create the lfn base directory
    lfn_base = law.wlcg.WLCGDirectoryTarget(
        dataset_key,
        fs=f"local",
    )
    # loop though files and interpret paths as lfns
    paths = [lfn_base.child(basename, type="f").path for basename in lfn_base.listdir(pattern="*.root")]

    return paths


def getGenTauDecayMode(prod: ak.Array):
    pids = prod.pdgId

    is_ele  = np.abs(pids) == 11
    is_muon = np.abs(pids) == 13
    is_charged = ((np.abs(pids) == 211) | (np.abs(pids) == 321))
    is_neutral = ((pids == 111) | (pids == 311) | (pids == 130) | (pids == 310))

    edecay = ak.sum(is_ele,  axis=-1) > 0
    mdecay = ak.sum(is_muon, axis=-1) > 0
    hdecay = (ak.sum(is_charged, axis=-1) > 0) | (ak.sum(is_neutral, axis=-1) >= 0)

    Nc = ak.sum(is_charged, axis=-1)
    Np = ak.sum(is_neutral, axis=-1)

    dm = ak.where(edecay, 
                  -1, 
                  ak.where(mdecay, 
                           -2, 
                           ak.where(hdecay, 
                                    (5 * (Nc - 1) + Np),
                                    -9)
                       )
              )

    return dm


def has_pt_greater_equal(obj_1, obj_2, pt_offset):
    """
    Checks if the pt of any object in vector1 is greater than or equal to 
    the pt of any object in vector2 plus an offset.
    
    Parameters:
        vector1 (awkward.Array): The pt values of the first vector (nested).
        vector2 (awkward.Array): The pt values of the second vector (nested).
        offset (float): The value to add to the pt of objects in vector2.
    
    Returns:
        awkward.Array: A boolean array with True if at least one object in vector1 
                       has pt greater than or equal to pt in vector2 + offset, otherwise False.
    """
    # Cartesian product to get all pairs
    pairs = ak.cartesian([obj_1, obj_2], axis=-1, nested=True)
    obj_1_br, obj_2_br = ak.unzip(pairs)
    return obj_1_br.pt >= (obj_2_br.pt + pt_offset)


@producer(
    uses={'event'},
    produces={
        "triggerID_*", "all_triggers_id",
    },
    exposed=False,
)
def hlt_path_matching(self: Producer, events: ak.Array, triggers: ak.Array, pair_objects: ak.Array, **kwags):
    # Initialize masks and dictionaries
    false_mask = ak.zeros_like(ak.local_index(events.event), dtype=np.bool_)
    single_electron_triggered = false_mask
    cross_electron_triggered  = false_mask
    single_muon_triggered = false_mask
    cross_muon_triggered  = false_mask
    cross_tau_triggered  = false_mask
    
    hlt_path_fired_e     = {}
    hlt_path_fired_etau  = {}
    hlt_path_fired_mu    = {}
    hlt_path_fired_mutau = {}
    hlt_path_fired_emu  = {}
    hlt_path_fired_tau   = {}
    
    matched_masks = {}
    trigger_ID = {}

    for trigger, trigger_fired, leg_masks in triggers.x.trigger_data:
        if trigger.has_tag("single_e"):
            single_ele_id = trigger.id
        if trigger.has_tag("single_mu"):
            single_mu_id = trigger.id
            
    # Perform each lepton election step separately per trigger
    for trigger, trigger_fired, leg_masks in triggers.x.trigger_data:
        
        is_single_el = trigger.has_tag("single_e")
        is_cross_el  = trigger.has_tag("cross_e_tau")
        is_single_mu = trigger.has_tag("single_mu")
        is_cross_mu  = trigger.has_tag("cross_mu_tau")
        is_cross_emu  = trigger.has_tag("cross_e_mu")
        is_cross_tau = trigger.has_tag("cross_tau_tau")
        
        dataset_name_tag = self.dataset_inst.keys[0]
          
        if (is_single_mu or is_cross_mu) & ('mutau' in pair_objects.keys()):
            muons = pair_objects.mutau.lep0
            taus = pair_objects.mutau.lep1
            
            if is_single_mu:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 13
                
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (muons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(muons.eta) < trigger.legs[0].min_eta)
                # muons to match with the trigger
                muons = muons[pt_mask0 & eta_mask0]
                
                # Delta R matching
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                single_mu_matches_leg0 = dr_matching
                single_mu_matches_leg0 = ak.any(ak.flatten(single_mu_matches_leg0, axis=-1), axis=1)
                single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                hlt_path_fired_mu[trigger.hlt_field] = ak.where(single_mu_matches_leg0, trigger.id, -1)
            
            elif is_cross_mu:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 13
                assert abs(trigger.legs[1].pdg_id) == 15
                
                dr_matching_mu = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                cross_mu_matches_leg0 = dr_matching_mu
    
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (muons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(muons.eta) < trigger.legs[0].min_eta)
                # pt requirement on the offline object before the trigger matching (leg 1)
                pt_mask1 = (taus.pt > trigger.legs[1].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 1)                 
                eta_mask1 = (abs(taus.eta) < trigger.legs[1].min_eta)
                
                # muons to match with the trigger
                muons = muons[pt_mask0 & eta_mask0]
                # taus to match with the trigger
                taus  = taus[pt_mask1 & eta_mask1]

                dr_matching_mu = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                cross_mu_matches_leg0 = dr_matching_mu
            
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_mu_tau_matches_leg1 = dr_matching_tau
                cross_mu_tau_matched = (ak.any(ak.flatten(cross_mu_matches_leg0, axis=-1), axis=1) & 
                                        ak.any(ak.flatten(cross_mu_tau_matches_leg1, axis=-1), axis=1))
                cross_muon_triggered = ak.where(trigger_fired & is_cross_mu, True, cross_muon_triggered)
                hlt_path_fired_mutau[trigger.hlt_field] = ak.where(cross_mu_tau_matched, trigger.id, -1)
        
            triggerID_mu    = hlt_path_fired(events, hlt_path_fired_mu)
            triggerID_mutau = hlt_path_fired(events, hlt_path_fired_mutau)
            matched_masks['mutau'] = ((triggerID_mu > 0) | (triggerID_mutau > 0))
            trigger_ID["triggerID_mu"]    = triggerID_mu
            trigger_ID["triggerID_mutau"] = triggerID_mutau
                

        if (is_single_el or is_cross_el) & ('etau' in pair_objects.keys()) :
            electrons =  pair_objects.etau.lep0
            taus =  pair_objects.etau.lep1
            if is_single_el:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 11
                
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (electrons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(electrons.eta) < trigger.legs[0].min_eta)
                
                # electrons to match with the trigger
                electrons = electrons[pt_mask0 & eta_mask0]           

                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                single_e_matches_leg0 =  dr_matching_e
                single_e_matches_leg0 = ak.any(ak.flatten(single_e_matches_leg0, axis=-1), axis=1)
                single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                hlt_path_fired_e[trigger.hlt_field] = ak.where(single_e_matches_leg0, trigger.id, -1)

            elif is_cross_el:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 11
                assert abs(trigger.legs[1].pdg_id) == 15
                
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (electrons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(electrons.eta) < trigger.legs[0].min_eta)
                # pt requirement on the offline object before the trigger matching (leg 1)
                pt_mask1 = (taus.pt > trigger.legs[1].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 1)                 
                eta_mask1 = (abs(taus.eta) < trigger.legs[1].min_eta) 
                
                # muons to match with the trigger
                electrons = electrons[pt_mask0 & eta_mask0]
                # taus to match with the trigger
                taus  = taus[pt_mask1 & eta_mask1]
                
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                cross_e_matches_leg0 = dr_matching_e
                
                dr_matching_tau = trigger_object_matching(taus, events.TrigObj[leg_masks[1]])
                cross_e_tau_matches_leg1 = dr_matching_tau
                
                cross_e_tau_matched = (ak.any(ak.flatten(cross_e_matches_leg0, axis=-1), axis=1) & 
                                       ak.any(ak.flatten(cross_e_tau_matches_leg1, axis=-1), axis=1))
                cross_electron_triggered = ak.where(trigger_fired & is_cross_el, True, cross_electron_triggered)
                hlt_path_fired_etau[trigger.hlt_field] = ak.where(cross_e_tau_matched, trigger.id, -1)
            
            triggerID_e     = hlt_path_fired(events, hlt_path_fired_e)
            triggerID_etau  = hlt_path_fired(events, hlt_path_fired_etau)
            matched_masks['etau'] = ((triggerID_e > 0) | (triggerID_etau > 0))
            trigger_ID["triggerID_e"]    = triggerID_e
            trigger_ID["triggerID_etau"] = triggerID_etau

        if (is_single_el or is_single_mu or is_cross_emu) & ('emu' in pair_objects.keys()) :
            electrons =  pair_objects.emu.lep0
            muons     =  pair_objects.emu.lep1
            
            if is_single_el:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 11
                
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (electrons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(electrons.eta) < trigger.legs[0].min_eta)
                
                # electrons to match with the trigger
                electrons = electrons[pt_mask0 & eta_mask0]           

                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[0]])
                single_e_matches_leg0 =  dr_matching_e
                single_e_matches_leg0 = ak.any(ak.flatten(single_e_matches_leg0, axis=-1), axis=1)
                single_electron_triggered = ak.where(trigger_fired & is_single_el, True, single_electron_triggered)
                hlt_path_fired_e[trigger.hlt_field] = ak.where(single_e_matches_leg0, trigger.id, -1)

            elif is_single_mu:
                assert trigger.n_legs == len(leg_masks) == 1
                assert abs(trigger.legs[0].pdg_id) == 13
                
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (muons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(muons.eta) < trigger.legs[0].min_eta)
                # muons to match with the trigger
                muons = muons[pt_mask0 & eta_mask0]
                
                # Delta R matching
                dr_matching = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                single_mu_matches_leg0 = dr_matching
                single_mu_matches_leg0 = ak.any(ak.flatten(single_mu_matches_leg0, axis=-1), axis=1)
                single_muon_triggered = ak.where(trigger_fired & is_single_mu, True, single_muon_triggered)
                hlt_path_fired_mu[trigger.hlt_field] = ak.where(single_mu_matches_leg0, trigger.id, -1)
                
            elif is_cross_emu:
                assert trigger.n_legs == len(leg_masks) == 2
                assert abs(trigger.legs[0].pdg_id) == 13
                assert abs(trigger.legs[1].pdg_id) == 11
                
                # pt requirement on the offline object before the trigger matching (leg 0)
                pt_mask0 = (muons.pt > trigger.legs[0].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 0)
                eta_mask0 = (abs(muons.eta) < trigger.legs[0].min_eta)
                # pt requirement on the offline object before the trigger matching (leg 1)
                pt_mask1 = (electrons.pt > trigger.legs[1].min_pt)
                # eta requirement on the offline object before the trigger matching (leg 1)                 
                eta_mask1 = (abs(electrons.eta) < trigger.legs[1].min_eta) 
                
                # muons to match with the trigger
                muons = muons[pt_mask0 & eta_mask0]
                # electrons to match with the trigger
                electrons  = electrons[pt_mask1 & eta_mask1]
                
                dr_matching_mu = trigger_object_matching(muons, events.TrigObj[leg_masks[0]])
                cross_mu_matches_leg0 = dr_matching_mu
                
                dr_matching_e = trigger_object_matching(electrons, events.TrigObj[leg_masks[1]])
                cross_e_matches_leg1 = dr_matching_e
                
                cross_mu_e_matched = (ak.any(ak.flatten(cross_mu_matches_leg0, axis=-1), axis=1) & 
                                       ak.any(ak.flatten(cross_e_matches_leg1, axis=-1), axis=1))
                cross_emu_triggered = ak.where(trigger_fired & is_cross_emu, True, cross_electron_triggered)
                hlt_path_fired_emu[trigger.hlt_field] = ak.where(cross_mu_e_matched, trigger.id, -1)          
            
            triggerID_e  = hlt_path_fired(events, hlt_path_fired_e)
            triggerID_mu = hlt_path_fired(events, hlt_path_fired_mu)
            triggerID_emu  = hlt_path_fired(events, hlt_path_fired_emu)

            
            # Trigger selection logic for emu channel:
            # -----------------------------------------
            # For Monte Carlo (MC) datasets:
            #   - Accept events if either of the triggers are fired:
            #       * HLT_IsoMu24 (represented by triggerID_mu > 0)
            #       * HLT_Ele30_WPTight_gsf (represented by triggerID_e > 0)
            #   This corresponds to:
            #       MC: HLT_IsoMu24 or HLT_Ele30_WPTight_gsf
            #
            # For real data, check if the 'emu' pair-object exists:
            #   - If the dataset is from the SingleMuon (or Muon_Run) channel:
            #       * Use HLT_IsoMu24 only (i.e., triggerID_mu > 0)
            #       This corresponds to:
            #           SingleMuon (Muon) : HLT_IsoMu24
            #
            #   - If the dataset is EGamma (NOT Muon_EGamma):
            #       * Accept events only if HLT_Ele30_WPTight_gsf is fired,
            #         and HLT_IsoMu24 is not fired.
            #         This is achieved by applying a mask that selects events where
            #         the electron trigger is active and the muon trigger is not.
            #       This corresponds to:
            #           EGamma : HLT_Ele30_WPTight_gsf and not HLT_IsoMu24
            # -----------------------------------------
            # For MC, set the default 'emu' mask to include events that fired any of the two triggers.
            if self.dataset_inst.is_mc:
                matched_masks['emu'] = ((triggerID_e > 0) | (triggerID_mu > 0))
            # For real data, apply the appropriate trigger selection if the 'emu' object is present.
            if 'emu' in pair_objects:
                # For SingleMuon or Muon_Run datasets, select events where the muon trigger (HLT_IsoMu24) fired.
                if ('SingleMuon' in dataset_name_tag) or ('Muon_Run' in dataset_name_tag):
                    matched_masks['emu'] = (triggerID_mu > 0)
                # For EGamma (MuonEG) datasets, select events where:
                # - The electron trigger (HLT_Ele30_WPTight_gsf) fired.
                # - The muon trigger (HLT_IsoMu24) did NOT fire.
                elif ('EGamma' in dataset_name_tag): #or ('MuonEG' in dataset_name_tag):
                    mask_e_not_mu_events = (triggerID_e == single_ele_id) & ~(triggerID_mu == single_mu_id)
                    # For events where the electron trigger is not isolated (i.e., doesn't meet the criteria), set its ID to -1.
                    triggerID_e = ak.where(mask_e_not_mu_events, triggerID_e, -1)
                    matched_masks['emu'] = (triggerID_e > 0)
                    
            trigger_ID["triggerID_e"]  = triggerID_e
            trigger_ID["triggerID_mu"] = triggerID_mu
            trigger_ID["triggerID_emu"] = triggerID_emu
            
        if is_cross_tau & ('tautau' in pair_objects.keys()) :
            assert trigger.n_legs == len(leg_masks) == 2
            assert abs(trigger.legs[0].pdg_id) == 15
            assert abs(trigger.legs[1].pdg_id) == 15
            
            taus1 =  pair_objects.tautau.lep0
            taus2 =  pair_objects.tautau.lep1
            
            # pt requirement on the offline object before the trigger matching (leg 0)
            pt_mask0 = (taus1.pt > trigger.legs[0].min_pt)
            # eta requirement on the offline object before the trigger matching (leg 0)
            eta_mask0 = (abs(taus1.eta) < trigger.legs[0].min_eta)
            # pt requirement on the offline object before the trigger matching (leg 1)
            pt_mask1 = (taus2.pt > trigger.legs[1].min_pt)
            # eta requirement on the offline object before the trigger matching (leg 1)                 
            eta_mask1 = (abs(taus2.eta) < trigger.legs[1].min_eta) 
            
            # muons to match with the trigger
            taus1 = taus1[pt_mask0 & eta_mask0]
            # taus to match with the trigger
            taus2 = taus2[pt_mask1 & eta_mask1]
            
            dr_matching_tau1 = trigger_object_matching(taus1, events.TrigObj[leg_masks[0]])
            dr_matching_tau2 = trigger_object_matching(taus2, events.TrigObj[leg_masks[1]])
            
            cross_tau1_matches_leg0 =  dr_matching_tau1
            cross_tau2_matches_leg1 =  dr_matching_tau2
            cross_tau_tau_matched = (ak.any(ak.flatten(cross_tau1_matches_leg0, axis=-1), axis=1) & 
                                     ak.any(ak.flatten(cross_tau2_matches_leg1, axis=-1), axis=1))
            cross_tau_triggered = ak.where(trigger_fired & is_cross_tau, True, cross_tau_triggered)
            hlt_path_fired_tau[trigger.hlt_field] = ak.where(cross_tau_tau_matched, trigger.id, -1)
        
            triggerID_tau   = hlt_path_fired(events, hlt_path_fired_tau)
            matched_masks['tautau'] = (triggerID_tau > 0)
            trigger_ID["triggerID_tau"]    = triggerID_tau  
    for column_name, array in trigger_ID.items():
        events = set_ak_column(events, column_name, array)
    matched_trigger_array = hlt_path_fired(events, trigger_ID)
    
    events = set_ak_column(events, "all_triggers_id", matched_trigger_array)

    return events, matched_masks


def find_fields_with_nan(awk_array):
    """
    Identifies fields in an Awkward Array that contain NaN values.

    Parameters:
        awk_array (ak.Array): Input Awkward Array.

    Returns:
        list: A list of field names that have NaN values.
    """
    fields_with_nan = []
    
    for field in awk_array.fields:
        field_data = awk_array[field]
        
        if ak.any(np.isnan(field_data)):
            fields_with_nan.append(field)
    
    return fields_with_nan