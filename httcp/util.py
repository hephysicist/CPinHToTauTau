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

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")


@deferred_column
def IF_NANO_V9(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version == 9 else None


@deferred_column
def IF_NANO_V11(self, func: ArrayFunction) -> Any | set[Any]:
    return self.get() if func.config_inst.campaign.x.version >= 10 else None


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
    any_match = ak.any(dr < threshold, axis=axis)

    return any_match


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


def enforce_hcand_type(hcand_pair_concat, field_type_dict):
    temp = {}
    for field, typename in field_type_dict.items():
        temp[field] = ak.enforce_type(ak.values_astype(hcand_pair_concat[field], typename), f"var * var * {typename}")
    hcand_array = ak.zip(temp)
    return hcand_array
    
