# coding: utf-8

"""
Selector related to additional MET filters.
"""

from __future__ import annotations

from columnflow.types import Iterable
from columnflow.selection import Selector, selector, SelectionResult
from columnflow.util import maybe_import
from columnflow.columnar_util import Route

ak = maybe_import("awkward")
np = maybe_import("numpy")


@selector(
    uses={"event","run","Flag*","Jet.{pt,eta,phi,mass,jetId,chEmEF,neEmEF}","PuppiMET.{pt,phi}"},
    
)
def met_filters_aux( 
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    mask = ak.ones_like(events.event, dtype=np.bool)
    if not self.dataset_inst.is_mc:
        run_mask = (events.run >= 362433) & (events.run <= 367144)
        BadCalibrationFilter_perjet_mask = (
            (events.PuppiMET.pt > 100)
            & (events.Jet.pt > 50)
            & (events.Jet.eta > -0.5)
            & (events.Jet.eta < -0.1)
            & (events.Jet.phi > -2.1)
            & (events.Jet.phi < -1.8)
            & ((events.Jet.neEmEF > 0.9) | (events.Jet.chEmEF > 0.9))
            & (abs(events.PuppiMET.delta_phi(events.Jet)) > 2.9)
        )
        BadCalibrationFilter_mask = ak.any(BadCalibrationFilter_perjet_mask, axis=1)
        mask = ak.where(run_mask, ~BadCalibrationFilter_mask, mask)

    return events, SelectionResult(steps={"bad_calib_filter": mask})
