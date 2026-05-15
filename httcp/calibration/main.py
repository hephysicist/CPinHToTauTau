"""
Calibration methods.
"""
import functools

from columnflow.calibration import Calibrator, calibrator
from httcp.calibration.jets import jec
from columnflow.calibration.cms.met import met_phi
from httcp.calibration.tau import tau_energy_scale
from httcp.calibration.electron import electron_smearing_scaling
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

import law

logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@calibrator(
    uses={
        jec, tau_energy_scale, deterministic_seeds, "PuppiMET.pt","Electron.phi", "Tau.phi", "Tau.pt"
    },
    produces={
        jec, tau_energy_scale, deterministic_seeds,"nJet", 
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:

    events = self[deterministic_seeds](events, **kwargs)
    
    # non_finite_mask = ~np.isfinite(events.PuppiMET.pt)
    #Jets variables before applying energy corrections
    events = set_ak_column_f32(events, "nJet", events.nJet)

    print("Performing Jet Energy Correction...")
    events = self[jec](events, **kwargs)
    
    #events = self[met_phi](events, **kwargs)
    #events = self[jer](events, **kwargs)
    #events = self[jets](events, **kwargs)
    
    if self.dataset_inst.is_mc: 
    #Apply tau energy scale correction
        print("Performing tau energy scale correction...")
        
        events = self[tau_energy_scale](events, do_syst=True, **kwargs)

    return events
