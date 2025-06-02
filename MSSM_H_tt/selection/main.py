# coding: utf-8

"""
Exemplary selection methods.
"""

from typing import Optional
from operator import and_
from functools import reduce

from collections import defaultdict, OrderedDict

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_stats
from columnflow.selection.cms.json_filter import json_filter
from columnflow.selection.cms.met_filters import met_filters

from columnflow.production.processes import process_ids
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from MSSM_H_tt.selection.physics_objects import muon_selection, electron_selection, tau_selection
from MSSM_H_tt.selection.trigger import trigger_selection
from MSSM_H_tt.selection.lepton_pair import pair_selection
from MSSM_H_tt.selection.lepton_veto import single_lepton_veto, second_lepton_veto, OC_lepton_veto
from MSSM_H_tt.selection.higgscand import new_higgscand, mask_nans
from MSSM_H_tt.selection.met_nanoAOD_filters import met_nanoAOD_filters
from MSSM_H_tt.production.aux_columns import channel_id
from MSSM_H_tt.selection.jets import jet_veto_map
from MSSM_H_tt.production.btag import btag_weight
from MSSM_H_tt.production.aux_columns import jets_taggable

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")


# exposed selectors
# (those that can be invoked from the command line)
@selector(
    uses={
        "event",
        # selectors / producers called within _this_ selector
        attach_coffea_behavior,
        json_filter,
        met_filters,
        mc_weight,
        process_ids,
        trigger_selection,
        muon_selection,
        electron_selection,
        tau_selection,
        pair_selection,
        channel_id,
        single_lepton_veto,
        second_lepton_veto,
        OC_lepton_veto,
        increment_stats,
        new_higgscand,
        mask_nans,
        jet_veto_map,
        btag_weight,
        jets_taggable,
        met_nanoAOD_filters,
    },
    produces={
        # selectors / producers whose newly created columns should be kept
        attach_coffea_behavior,
        json_filter,
        met_filters,
        mc_weight,
        process_ids,
        trigger_selection,
        muon_selection,
        electron_selection,
        tau_selection,
        pair_selection,
        channel_id,
        single_lepton_veto,
        second_lepton_veto,
        OC_lepton_veto,
        increment_stats,
        new_higgscand,
        mask_nans,
        jet_veto_map,
        btag_weight,
        jets_taggable,
        met_nanoAOD_filters,
        "category_ids",
        "OC_lepton_veto",
    },
    exposed=True,
)
def main(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    # ensure coffea behaviors are loaded
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # filter bad data events according to golden lumi mask
    if self.dataset_inst.is_data:
        events, json_filter_results = self[json_filter](events, **kwargs)
        results += json_filter_results

    # trigger selection
    events, trigger_results = self[trigger_selection](events, **kwargs)
    results += trigger_results

    # met filter selection
    events, met_filter_results = self[met_filters](events, **kwargs)
    results += met_filter_results

    # muon selection
    # e.g. mu_idx: [ [0,1], [], [1], [0], [] ]
    events, muon_results, good_muon_indices, mu_emu_indices, single_veto_muon_indices, second_mu_veto_emu_indices, OC_veto_muon_indices = self[muon_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)
    
    # electron selection
    # e.g. ele_idx: [ [], [0,1], [], [], [1,2] ]
    events, ele_results, good_electron_indices, ele_emu_indices, single_veto_electron_indices, second_ele_veto_emu_indices, OC_veto_electron_indices = self[electron_selection](events,
                                                                                                           call_force=True,
                                                                                                           **kwargs)

    # tau selection
    # e.g. tau_idx: [ [1], [0,1], [1,2], [], [0,1] ]
    events, good_tau_indices = self[tau_selection](events,
                                                    call_force=True,
                                                    **kwargs)
    
    mutau_indices_pair = self[pair_selection](events,
                                              'mutau',
                                              good_muon_indices,
                                              good_tau_indices)
    etau_indices_pair = self[pair_selection](events,
                                             'etau',
                                             good_electron_indices,
                                             good_tau_indices)
    emu_indices_pair = self[pair_selection](events,
                                             'emu',
                                             ele_emu_indices,
                                             mu_emu_indices)
    
    tautau_indices_pair = self[pair_selection](events,
                                              'tautau',
                                              good_tau_indices,
                                              good_tau_indices)
    
    pair_idxs = DotDict.wrap({
        'etau'  : etau_indices_pair,
        'mutau' : mutau_indices_pair,
        'emu'   : emu_indices_pair,
        'tautau': tautau_indices_pair
    })
    
    raw_dilepton_mask = ak.ones_like(events.event, dtype=np.bool_)
    for channel in self.config_inst.channels.names():
        raw_dilepton_mask = raw_dilepton_mask | eval(f'(ak.num({channel}_indices_pair.lep0, axis=1) > 0)')
    
    results += SelectionResult(
    steps = {
        "has_at_least_2_leptons" : raw_dilepton_mask,
    })
    events, hcand_res = self[new_higgscand](events,
                                 trigger_results,
                                 pair_idxs,
                                 domatch=True,
                                 **kwargs)
    results += hcand_res
    
    #produce channel id column (legacy)
    events = self[channel_id](events)

    # Single lepton veto
    # it is only applied on the events with one higgs candidate only
    events, single_lepton_veto_results = self[single_lepton_veto](events,
                                                                  single_veto_electron_indices,
                                                                  single_veto_muon_indices)
    results += single_lepton_veto_results

    # # Additional lepton veto
    if self.config_inst.channels.names()[0] == 'emu':
        events, second_lepton_veto_results = self[second_lepton_veto](events,
                                                                      second_ele_veto_emu_indices,
                                                                      second_mu_veto_emu_indices)
    else:
        events, second_lepton_veto_results = self[second_lepton_veto](events,
                                                                      single_veto_electron_indices,
                                                                      single_veto_muon_indices)
    results += second_lepton_veto_results

    # Opposite Charge (OC) lepton pair veto
    events, OC_lepton_veto_results = self[OC_lepton_veto](events,
                                                        OC_veto_electron_indices,
                                                        OC_veto_muon_indices)
    
    #Check arrays for np.nan values and mask them
    events, nan_mask_res = self[mask_nans](events)
    results += nan_mask_res
    
    # additional jet veto map, vetoing entire events
    if self.has_dep(jet_veto_map):
        events, jet_veto_map_result = self[jet_veto_map](events, **kwargs)
        results += jet_veto_map_result
    
    if self.dataset_inst.is_data:
        events, met_nanoAOD_filters_result = self[met_nanoAOD_filters](events, **kwargs)
        results += met_nanoAOD_filters_result
    
    # combined event selection after all steps
    event_sel = reduce(and_, results.steps.values())
    
    results.event = event_sel
    
    events = self[jets_taggable](events, **kwargs) 
    # add the mc weight
    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)
        events = self[btag_weight](events, do_syst=True, **kwargs)
        events_b_weight = self[mc_weight](events, **kwargs)
        w_event_btag = events.btag_weight_nom
    events = self[process_ids](events, **kwargs)
    events = set_ak_column(events, 'category_ids', ak.ones_like(events.event, dtype=np.uint8))
    
    njets_taggable = events.n_jets_tag
    
    weight_map = {
        "num_events": Ellipsis,
        "num_events_selected": event_sel,
    }
    group_map = {}
    if self.dataset_inst.is_mc:
        weight_map = {
            **weight_map,
            # mc weight for all events
            "sum_mc_weight": (events.mc_weight, Ellipsis),
            "sum_mc_weight_selected": (events.mc_weight, results.event),
            "sum_mc_weight_selected_with_btag_weight": (events.mc_weight*w_event_btag, results.event),
        }
        group_map = {
            # per process
            "process": {
                "values": events.process_id,
                "mask_fn": (lambda v: events.process_id == v),
            },
            # per channel
            "channel": {
                "values": events.channel_id,
                "mask_fn": (lambda v: events.channel_id == v),
            },
        
        }
        # per jet multiplicity
        if njets_taggable is not None:
            group_map["njets_taggable"] = {
                "values": njets_taggable,
                "mask_fn": (lambda v: njets_taggable == v),
            }
    events, results = self[increment_stats](
        events, results, stats, weight_map=weight_map, group_map=group_map, **kwargs)
    
    return events, results