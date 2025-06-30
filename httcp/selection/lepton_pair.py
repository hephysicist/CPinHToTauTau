# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""
import copy
from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from httcp.util import get_lep_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

def select_best_pair(
        lep0: ak.Array,
        lep1: ak.Array,
        vars4sorting: dict,
)->ak.Array:
    def leading_lep(lep):
        return ak.firsts(lep[:,:1], axis=1)
    empty_lep = ak.zeros_like(lep0.obj_idx, dtype=np.int32)[..., :0]
    lep0_idx = empty_lep
    lep1_idx = empty_lep
    for (var_str, sort_dir) in vars4sorting.items():
        the_var = eval(var_str)
        sorted_idx = ak.argsort(the_var, ascending=(sort_dir == 'ascending'))
        sorted_var = the_var[sorted_idx]
        vars_not_equal = ak.fill_none((ak.firsts(sorted_var[:,:1], axis=1) != ak.firsts(sorted_var[:,1:2], axis=1)),False)
        lep0_leading_idx = lep0.obj_idx[sorted_idx[:,:1]]
        lep0_idx = ak.where(vars_not_equal,lep0_leading_idx, lep0_idx)
        lep1_leading_idx = lep1.obj_idx[sorted_idx[:,:1]]
        lep1_idx = ak.where(vars_not_equal,lep1_leading_idx, lep1_idx)
    pairs = ak.zip({'lep0': lep0_idx,
                    'lep1': lep1_idx})
    return pairs

@selector(
    uses=
        {f"Electron.{var}" for var in [
            "pt", "eta", "phi", "mass","pfRelIso03_all", "rawIdx"]
        } | {
            f"Muon.{var}" for var in [
                "pt","eta","phi","mass", "pfRelIso04_all"] 
        } | {
            f"Tau.{var}" for var in [
                "pt","eta","phi","mass","charge", 
                "idDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSe","idDeepTau2018v2p5VSmu","rawDeepTau2018v2p5VSjet"
            ] 
        },
    exposed=False,
)
def pair_selection(
        self: Selector,
        events: ak.Array,
        channel: str,
        lep0_idxs: ak.Array,
        lep1_idxs: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:
    print(f'Selecting dilepton pairs for {channel}')
    ch_objects = self.config_inst.x.ch_objects
    presel_leps = []
    for idx in range(2):
        the_lep = eval(f"events[ch_objects[channel].lep{idx}][lep{idx}_idxs]")
        the_lep['obj_idx'] = eval(f"lep{idx}_idxs")
        presel_leps.append(the_lep)
    pairs = ak.cartesian(presel_leps, axis=1)
    lep0, lep1 = ak.unzip(pairs)
    #Can be the case that this line is not needed, but it is better to explicitly define 4-vectors so the methods of coffea work properly
    [lep0_p4,lep1_p4] = [get_lep_p4(the_lep) for the_lep in [lep0,lep1]]
    if channel =='mutau': 
        vars4sorting = {'lep0.pfRelIso04_all'           : 'ascending', 
                        'lep0.pt'                       : 'descending',
                        'lep1.rawDeepTau2018v2p5VSjet'  : 'ascending',
                        'lep1.pt'                       : 'descending'}
        pair_cuts = {"dr_0p5"        : lep0_p4.delta_r(lep1_p4) > 0.5, 
                     "invmass_40"    : (lep0_p4 + lep1_p4).mass > 40, }
    if channel =='etau': 
        vars4sorting = {'lep0.pfRelIso03_all'           : 'ascending', 
                        'lep0.pt'                       : 'descending',
                        'lep1.rawDeepTau2018v2p5VSjet'  : 'ascending',
                        'lep1.pt'                       : 'descending'}
        pair_cuts = {"dr_0p5"        : lep0_p4.delta_r(lep1_p4) > 0.5}
    if channel =='tautau': 
        vars4sorting = {'lep0.rawDeepTau2018v2p5VSjet'  : 'ascending', 
                        'lep1.rawDeepTau2018v2p5VSjet'  : 'ascending',
                        'lep0.pt'                       : 'descending',
                        'lep1.pt'                       : 'descending'}
        pair_cuts = {"is_pt_40"   : (lep0.pt > 40) & (lep1.pt > 40),
                     "eta_2p1"    : (np.abs(lep0.eta) < 2.1) & (np.abs(lep1.eta) < 2.1),
                     "dr_0p5"     : lep0.delta_r(lep1) > 0.5}
    mask = ak.ones_like(lep0_p4.pt, dtype=np.bool_)
    for cut in pair_cuts.values():
        mask = mask & cut
    [presel_lep0,presel_lep1] = [the_lep[mask] for the_lep in [lep0,lep1]]
    has_multiple_pairs = ak.num(lep0.pt, axis=1) > 1
    pair_idxs = ak.where(has_multiple_pairs,
                         select_best_pair(presel_lep0,presel_lep1,vars4sorting),
                         ak.zip({'lep0': presel_lep0.obj_idx,
                                 'lep1': presel_lep1.obj_idx}))
    return pair_idxs
