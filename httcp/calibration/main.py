# coding: utf-8

"""
main calibration script
"""

from columnflow.calibration import Calibrator, calibrator
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from httcp.calibration.tau import tau_energy_scale

np = maybe_import("numpy")
ak = maybe_import("awkward")


@calibrator(
    uses={
        deterministic_seeds,
        tau_energy_scale,
    },
    produces={
        deterministic_seeds,
        tau_energy_scale,
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    events = self[deterministic_seeds](events, **kwargs)
    if self.dataset_inst.is_mc: 
    #Apply tau energy scale correction
        events = self[tau_energy_scale](events, **kwargs)
    return events
