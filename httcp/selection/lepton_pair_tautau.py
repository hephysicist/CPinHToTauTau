# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column


np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



def get_sorted_pair(
        dtrpairs: ak.Array,
        dtrpairindices: ak.Array,
)->ak.Array:
    # redundant, because taus were sorted by the deeptau before
    sorted_idx = ak.argsort(dtrpairs["0"].rawDeepTau2018v2p5VSjet, ascending=False)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # if the deep tau val of tau-0 is the same for the first two pair
    where_same_iso_1 = (
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # if so, sort the pairs according to the deep tau of the 2nd tau
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["1"].rawDeepTau2018v2p5VSjet, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # if the deep tau val of tau-1 is the same for the first two pair 
    where_same_iso_2 = (
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["1"].rawDeepTau2018v2p5VSjet[:,1:2], axis=1)
    )

    # sort them with the pt of the 1st tau
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # check if the first two pairs have the second tau with same rawDeepTau2017v2p1VSjet
    where_same_pt_1 = (
        ak.firsts(dtrpairs["0"].pt[:,:1], axis=1)
        ==
        ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1)
    )

    # if so, sort the taus with their pt
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)

    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    lep1idx = ak.singletons(ak.firsts(dtrpairindices["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(dtrpairindices["1"], axis=1))
    # print(f"lep1 pt: {lep1.pt}")
    # print(f"lep2 pt: {lep2.pt}")
    dtrpair    = ak.concatenate([lep1, lep2], axis=1) 
    dtrpairidx = ak.concatenate([lep1idx, lep2idx], axis=1)     

    return dtrpairidx



@selector(
    uses={
        "Tau.pt", "Tau.eta", "Tau.phi", "Tau.mass",
        "Tau.charge", "Tau.rawDeepTau2018v2p5VSjet",
    },
    exposed=False,
)
def tautau_selection(
        self: Selector,
        events: ak.Array,
        lep_indices: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:

    #from IPython import embed; embed()

    # Sorting leps [Tau] by deeptau [descending]
    lep_sort_key       = events.Tau[lep_indices].rawDeepTau2018v2p5VSjet
    lep_sorted_indices = ak.argsort(lep_sort_key, axis=-1, ascending=False)
    lep_indices        = lep_indices[lep_sorted_indices]

    leps_pair        = ak.combinations(events.Tau[lep_indices], 2, axis=1)
    lep_indices_pair = ak.combinations(lep_indices, 2, axis=1)
    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)
    lep1_idx, lep2_idx = ak.unzip(lep_indices_pair)

    preselection = {
        "tautau_is_pt_40"      : (lep1.pt > 40) & (lep2.pt > 40),
        "tautau_is_eta_2p1"    : (np.abs(lep1.eta) < 2.1) & (np.abs(lep2.eta) < 2.1),
        "tautau_is_os"         : (lep1.charge * lep2.charge) < 0,
        "tautau_dr_0p5"        : (1*lep1).delta_r(1*lep2) > 0.5,  #deltaR(lep1, lep2) > 0.5,
    }

    good_pair_mask = lep1_idx >= 0
    pair_selection_steps = {}
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = good_pair_mask
        
    leps_pair_sel = leps_pair[good_pair_mask]
    lep_indices_pair_sel = lep_indices_pair[good_pair_mask]

    lep1idx = ak.singletons(ak.firsts(lep_indices_pair_sel["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(lep_indices_pair_sel["1"], axis=1))

    lep_indices_pair_sel_single = ak.concatenate([lep1idx, lep2idx], axis=1)

    where_many   = ak.num(lep_indices_pair_sel, axis=1) > 1
    pair_indices = ak.where(where_many, 
                            get_sorted_pair(leps_pair_sel,
                                            lep_indices_pair_sel),
                            lep_indices_pair_sel_single)

    return SelectionResult(
        aux = pair_selection_steps,
    ), pair_indices
