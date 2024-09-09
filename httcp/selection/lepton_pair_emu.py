# coding: utf-8

"""
Prepare h-Candidate from SelectionResult: selected lepton indices & channel_id [trigger matched] 
"""

from typing import Optional
from columnflow.selection import SelectionResult
from columnflow.selection.util import create_collections_from_masks
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column

from hcp.util import invariant_mass, deltaR, transverse_mass

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


def get_sorted_pair(
        dtrpairs: ak.Array,
        dtrpairindices: ak.Array
)->ak.Array:

    sorted_idx = ak.argsort(dtrpairs["0"].pfRelIso03_all, ascending=True)
    # Sort the pairs based on pfRelIso03_all of the first object in each pair
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # Extract the pfRelIso03_all values for the first object in each pair
    lep1_pfRelIso03_all = dtrpairs["0"].pfRelIso03_all
    # Check if the pfRelIso03_all values are the same for the first two objects in each pair
    where_same_iso_1 = ak.fill_none(
        (
            ak.firsts(dtrpairs["0"].pfRelIso03_all[:,:1], axis=1) 
            ==
            ak.firsts(dtrpairs["0"].pfRelIso03_all[:,1:2], axis=1)
        ), False
    )
    
    # Sort the pairs based on pt if pfRelIso03_all is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_1,
                          ak.argsort(dtrpairs["0"].pt, ascending=False),
                          sorted_idx)

    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # Check if the pt values are the same for the first two objects in each pair    
    where_same_pt_1 = ak.fill_none(
        (
            ak.firsts(dtrpairs["0"].pt[:,:1], axis=1)
            ==
            ak.firsts(dtrpairs["0"].pt[:,1:2], axis=1)
        ), False
    )
    
    # if so, sort the pairs with tau rawDeepTau2017v2p1VSjet
    sorted_idx = ak.where(where_same_pt_1,
                          ak.argsort(dtrpairs["1"].pfRelIso03_all, ascending=False),
                          sorted_idx)
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]
    # check if the first two pairs have taus with same rawDeepTau2017v2p1VSjet
    where_same_iso_2 = ak.fill_none(
        (
            ak.firsts(dtrpairs["1"].pfRelIso03_all[:,:1], axis=1)
            ==
            ak.firsts(dtrpairs["1"].pfRelIso03_all[:,1:2], axis=1)
        ), False
    )
    
    # Sort the pairs based on pt if rawDeepTau2017v2p1VSjet is the same for the first two objects
    sorted_idx = ak.where(where_same_iso_2,
                          ak.argsort(dtrpairs["1"].pt, ascending=False),
                          sorted_idx)
    # finally, the pairs are sorted
    dtrpairs = dtrpairs[sorted_idx]
    dtrpairindices = dtrpairindices[sorted_idx]

    # Extract the first object in each pair (lep1) and the second object (lep2)
    lep1 = ak.singletons(ak.firsts(dtrpairs["0"], axis=1))
    lep2 = ak.singletons(ak.firsts(dtrpairs["1"], axis=1))

    lep1idx = ak.singletons(ak.firsts(dtrpairindices["0"], axis=1))
    lep2idx = ak.singletons(ak.firsts(dtrpairindices["1"], axis=1))

    # Concatenate lep1 and lep2 to create the final dtrpair
    dtrpair    = ak.concatenate([lep1, lep2], axis=1)
    dtrpairidx = ak.concatenate([lep1idx, lep2idx], axis=1)

    return dtrpairidx



@selector(
    uses={
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
        "Electron.charge", "Electron.pfRelIso03_all",
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass",
        "Muon.charge", "Muon.pfRelIso03_all",
        "PuppiMET.pt", "PuppiMET.phi",
    },
    exposed=False,
)
def emu_selection(
        self: Selector,
        events: ak.Array,
        lep1_indices: ak.Array,
        lep2_indices: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult, ak.Array]:

    lep_indices_pair = ak.cartesian([lep1_indices, 
                                     lep2_indices], axis=1)
    leps_pair        = ak.cartesian([events.Electron[lep1_indices], 
                                     events.Tau[lep2_indices]], axis=1)
    
    # pair of leptons: probable higgs candidate -> leps_pair
    # and their indices                         -> lep_indices_pair 
    lep1, lep2 = ak.unzip(leps_pair)
    lep1_idx, lep2_idx = ak.unzip(lep_indices_pair)

    preselection = {
        "is_os"         : (lep1.charge * lep2.charge) < 0,
        "dr_0p5"        : deltaR(lep1, lep2) > 0.5,
        #"mT_50"         : transverse_mass(lep1, events.PuppiMET) < 50
    }

    good_pair_mask = lep1_idx >= 0
    pair_selection_steps = {}
    for cut in preselection.keys():
        good_pair_mask = good_pair_mask & preselection[cut]
        pair_selection_steps[cut] = ak.sum(preselection[cut], axis=1) > 0
        

    leps_pair_sel = leps_pair[good_pair_mask]
    lep_indices_pair_sel = lep_indices_pair[good_pair_mask]

    where_many   = ak.num(lep_indices_pair_sel, axis=1) > 1
    pair_indices = ak.where(where_many, 
                            get_sorted_pair(leps_pair_sel,
                                            lep_indices_pair_sel),
                            lep_indices_pair_sel)


    return events, SelectionResult(
        steps = pair_selection_steps,
    ), pair_indices
