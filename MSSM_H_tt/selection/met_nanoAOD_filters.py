# coding: utf-8
from __future__ import annotations

import law
import math

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view, optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")


logger = law.logger.get_logger(__name__)

# The NanoAOD does not have enough info to rerun the filter. So please apply the following recipe:

# Reject the event if PuppiMET _pt > 100 GeV and there is at least one jet (AK4) which has

#         pT > 50 GeV,
#         η within -0.5 to -0.1,
#         Φ within -2.1 to -1.8,
#         Neutral EM energy fraction or charged EM energy fraction (branch names: Jet_neEmEF, Jet_chEmEF) > 0.9
#         ΔΦ(PuppiMET _phi, jet) > 2.9 
#     Apply it only for RunNumbers in the range 362433 to 367144 which belong to later part of 2022 and early 2023.
#     DO NOT apply jet ID (branch: Jet_jetId) on the jets while implementing this recipe.
#     The effect of this recipe on good events is very small (<0.2%) and it is not simulated in MC. So, the recipe is not recommended for MC.
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
@selector(
    uses={
        "Jet.{chEmEF,neEmEF}", 
        "PuppiMET.{pt,phi}",
    },
)
def met_nanoAOD_filters(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:

    jet = events.Jet
    met = events.PuppiMET
    jet_mask = (
        (jet.pt > 50) &
        (jet.eta > -0.5)&
        (jet.eta < -0.1)&
        (jet.phi > -2.1)&
        (jet.phi < -1.8)&
        ((jet.chEmEF< 0.9)|(jet.neEmEF< 0.9))&
        ((met.phi - jet.phi) >= 2.9)
        )
    met_filter_mask = (events.run >= 362433) & (events.run <= 367144) & (met.pt >= 100) & ak.any(jet_mask,axis=1)

    # create the selection result
    results = SelectionResult(
        steps={"met_nanoAOD_filters": ~met_filter_mask},
    )

    return events, results
