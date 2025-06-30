# coding: utf-8

"""
Selection modules for jets.
"""

from __future__ import annotations

import law
import math

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view, optional_column as optional
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
np = maybe_import("numpy")
ak = maybe_import("awkward")


logger = law.logger.get_logger(__name__)


@selector(
    uses={
        "Jet.{pt,eta,phi,mass,jetId,chEmEF,neEmEF}", 
        "Muon.{pt,eta,phi,mass,isPFcand}",
        optional("Jet.puId"),
    },
    produces={"Jet.veto_map_mask"},
    get_veto_map_file=(lambda self, external_files: external_files.jet_veto_map),
)
def jet_veto_map(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    Selector that applies the Jet Veto Map to the jets and stores the result as a new column ``Jet.veto_maps``.
    Additionally, the ``jet_veto_map`` step is added to the SelectionResult that masks events containing
    jets from the veto map, which is the recommended way to use the veto map.
    For users that only want to remove the jets from the veto map, the ``veto_map_jets`` object
    is added to the SelectionResult.

    Requires an external file in the config
    under ``jet_veto_map``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "jet_veto_map": ("/afs/cern.ch/user/m/mfrahm/public/mirrors/jsonpog-integration-a332cfa/POG/JME/2022_Summer22EE/jetvetomaps.json.gz", "v1"),  # noqa
        })

    *get_veto_map_file* can be adapted in a subclass in case it is stored differently in the external files.

    documentation: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
    """
    jet = events.Jet
    muon = events.Muon[events.Muon.isPFcand]
    # loose jet selection
    jet_mask = (
        (jet.pt > 15) &
        (jet.jetId >= 2) &  # tight id 
        ((jet.chEmEF + jet.neEmEF) < 0.9) &
        ak.all(jet.metric_table(muon) >= 0.2, axis=2) &
        (np.abs(jet.eta) < 5.191)    # https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/JME_2022_Prompt_jetvetomaps.html
    )
    #Create 
    presel_jet = ak.mask(jet,jet_mask)
    
    veto_args = ["jetvetomap", presel_jet.eta, presel_jet.phi]
    veto_mask = ak.fill_none(self.veto_map.evaluate(*veto_args),False)
    jet_veto_per_evt = ~ak.any(veto_mask, axis=1)
    events = set_ak_column(events, "Jet.veto_map_mask", veto_mask)
    # create the selection result
    results = SelectionResult(
        steps={"jet_veto_map": jet_veto_per_evt},
    )
    return events, results


@jet_veto_map.requires
def jet_veto_map_requires(self: Selector, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@jet_veto_map.setup
def jet_veto_map_setup(
    self: Selector,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]

    # create the corrector
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_veto_map_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )
    keys = list(correction_set.keys())
    if len(keys) != 1:
        raise ValueError(f"Expected exactly one correction in the file, got {len(keys)}")

    self.veto_map = correction_set[keys[0]]
