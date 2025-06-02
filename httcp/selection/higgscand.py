# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched]
"""
import functools
from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.util import DotDict

from httcp.util import hlt_path_matching, hlt_path_fired

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@selector(
    uses={
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "rawDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSjet", "idDeepTau2018v2p5VSe", "idDeepTau2018v2p5VSmu",
                "rawIdx", "IPx", "IPy", "IPz", "ip_sig", "jetIdx"]
    } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge",
                "rawIdx","IPx", "IPy", "IPz","ip_sig", "jetIdx"
            ] 
    } | {
            f"Electron.{var}" for var in [
                "pt","eta","phi","mass","dxy","dz", "charge", 
                "pfRelIso03_all", "rawIdx", "IPx", "IPy", "IPz","ip_sig", "jetIdx"
            ] 
        } | {optional("Tau.genPartFlav")} | {hlt_path_matching},
    produces={
        'hcand_*'
    } | {hlt_path_matching},
    exposed=False,
)
def new_higgscand(
        self: Selector,
        events: ak.Array,
        trigger_results: SelectionResult,
        pair_idxs: dict,
        domatch: Optional[bool] = False,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    steps = {}
    pair_objects = {}
    n_pairs_prematch = ak.zeros_like(events.event)
    for the_ch in channels:
        obj_idx = pair_idxs[the_ch]
        pair_objects[the_ch] = {}
        for lep in obj_idx.fields : 
            pair_objects[the_ch][lep] = events[ch_objects[the_ch][lep]][obj_idx[lep]]
        n_pairs_prematch = n_pairs_prematch + ak.num(pair_objects[the_ch]['lep0'].pt>0, axis = 1)
    steps['selected_hcand'] = n_pairs_prematch > 0
    pair_objects = DotDict.wrap(pair_objects)
    
    if domatch: 
        n_pairs_postmatch = ak.zeros_like(events.event)
        events, matched_masks = self[hlt_path_matching](events, trigger_results, pair_objects, **kwargs)
        # Define arrays in a dictionary for easy management
        hcands = {}
        for the_ch in channels:
            matched_pairs = {}
            for lep_idx in range(2):
                the_lep = pair_objects[the_ch][f'lep{lep_idx}']
                matched_pairs[f'lep{lep_idx}'] = ak.where(matched_masks[the_ch], the_lep, the_lep[:,:0])
            events = set_ak_column(events, f"hcand_{the_ch}",  ak.zip(matched_pairs))
            n_pairs_postmatch =  n_pairs_postmatch + ak.num(events[f"hcand_{the_ch}"].lep0.pt>0, axis = 1)
        steps['selected_hcand_trigmatch'] = n_pairs_postmatch > 0
        steps['single_hcand'] = (n_pairs_postmatch == 1)
    else:
        for the_ch in channels:
            events = set_ak_column(events, f"hcand_{the_ch}",  ak.zip(pair_objects[the_ch]))
            steps['single_hcand'] = n_pairs_prematch == 1
    return events, SelectionResult(steps=steps)

@selector(
    uses={
        'hcand_*'
    },
    exposed=False,
)
def mask_nans(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
       
        hcand = events[f'hcand_{ch_str}']
        for lep_str in [field for field in hcand.fields if 'lep' in field]:
            lep = hcand[lep_str]
            for field in lep.fields:
                field_mask = (ak.any(np.isnan(lep[field]),axis=1))
                mask = mask | field_mask
                if ak.any(field_mask):
                    evt_index = events.event[field_mask]
                    print(f'Found nans:\n{ch_str}.{lep_str}.{field}:')
                    for evt in evt_index:
                        print(f'evt idx: {evt}')
    return events, SelectionResult(
        steps={
            'nans_removed': ~mask
        }
        )
    
