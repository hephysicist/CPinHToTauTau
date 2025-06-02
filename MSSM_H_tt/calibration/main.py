"""
Calibration methods.
"""
import functools

from columnflow.calibration import Calibrator, calibrator
from MSSM_H_tt.calibration.jets import jec
from MSSM_H_tt.calibration.met import met_phi
from MSSM_H_tt.calibration.tau import tau_energy_scale
from MSSM_H_tt.calibration.electron import electron_smearing_scaling
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
        jec, tau_energy_scale, electron_smearing_scaling, deterministic_seeds, "Electron.phi", "Tau.phi", "Tau.pt", "run", "luminosityBlock", "event",
    },
    produces={
        jec, tau_energy_scale, electron_smearing_scaling, deterministic_seeds, "Jet.pt_no_jec", "Jet.eta_no_jec",
        "Jet.phi_no_jec", "Jet.mass_no_jec", "PuppiMET.pt_no_jec", "PuppiMET.phi_no_jec", "nJet", "Electron.pt_no_scaling_smearing", 
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    
    events = self[deterministic_seeds](events, **kwargs)
    
    non_finite_mask = ~np.isfinite(events.PuppiMET.pt)
    
    #Jets variables before applying energy corrections
    events = set_ak_column_f32(events, "Jet.pt_no_jec", events.Jet.pt)
    events = set_ak_column_f32(events, "Jet.phi_no_jec", events.Jet.phi)
    events = set_ak_column_f32(events, "Jet.eta_no_jec", events.Jet.eta)
    events = set_ak_column_f32(events, "Jet.mass_no_jec", events.Jet.mass)
    events = set_ak_column_f32(events, "nJet", events.nJet)

    #PuppiMET variables before applying energy corrections
    events = set_ak_column_f32(events, "PuppiMET.pt_no_jec", events.PuppiMET.pt)
    events = set_ak_column_f32(events, "PuppiMET.phi_no_jec", events.PuppiMET.phi)

    print("Performing Jet Energy Correction...")
    events = self[jec](events, **kwargs)
    events = set_ak_column_f32(events, "Electron.pt_no_scaling_smearing", events.Electron.pt)
    
    if self.config_inst.x.year==2022: # scaling and smearing is not available for 2023
        print("Performing electron scaling and smearing correction...")
        events = self[electron_smearing_scaling](events, **kwargs)
        print("Electron scaling and smearing correction...SUCCEDED")
    #events = self[met_phi](events, **kwargs)
    #events = self[jer](events, **kwargs)
    #events = self[jets](events, **kwargs)
    if self.dataset_inst.is_mc & (self.config_inst.channels.names()[0]!= 'emu'): 
    #Apply tau energy scale correction
        print("Performing tau energy scale correction...")
        
        events = self[tau_energy_scale](events, **kwargs)
   
    return events
