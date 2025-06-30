# coding: utf-8

"""
Exemplary calibration methods.
"""
import functools
import itertools

from columnflow.calibration import Calibrator, calibrator
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import, InsertableDict
from columnflow.columnar_util import set_ak_column, flat_np_view


np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

###############################################################
### ELECTRON SMEARING ONLY FOR MC AND SCALING ONLY FOR DATA ###
###############################################################

@calibrator(
    uses={"Electron.pt","Electron.r9","Electron.eta", "Electron.seedGain", "Electron.scEtOverPt", "run", 
    },
    produces={
        "Electron.pt",
    },
    mc_only=False,
)
def electron_smearing_scaling(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    # fail when running on data
    if self.dataset_inst.is_data:    
        syst = "total_correction" 
        gain = events.Electron.seedGain #int
        run = events.run #real
        eta = events.Electron.eta #real
        r9 = events.Electron.r9 #real
        et = events.Electron.scEtOverPt #real

        electron_pt = events.Electron.pt
        
        #Create get energy scale correction for each tau
        electron_scaling_nom = np.ones_like(electron_pt, dtype=np.float32)

        if self.config_inst.x.year in [2022,2023]:
            electron_scaling_args = lambda events, syst: ( syst,
                                                            gain,
                                                            run,
                                                            eta,
                                                            r9,
                                                            et,
                                                            )

        electron_scaling_nom = self.electron_scaling_corrector.evaluate(*electron_scaling_args(events, syst))

        events = set_ak_column_f32(events, "Electron.pt", electron_pt * electron_scaling_nom)
    
    elif self.dataset_inst.is_mc:
        
        syst = "rho"  
        eta =events.Electron.eta 
        r9 = events.Electron.r9 
        electron_pt = events.Electron.pt
        #Create get energy scale correction for each electron
        electron_smearing_nom = np.ones_like(electron_pt, dtype=np.float32)

        if self.config_inst.x.year in [2022,2023]:
            electron_smearing_args = lambda events, syst: (syst,
                                                        eta,
                                                        r9)

        electron_smearing_nom = self.electron_smearing_corrector.evaluate(*electron_smearing_args(events, syst))
        rng = np.random.default_rng(seed=125)  
        
        electron_smearing_nom_flat      = np.asarray(ak.flatten(electron_smearing_nom))
        electron_pt      = np.asarray(ak.flatten(events.Electron.pt))
        arr_shape = ak.num(events.Electron.pt, axis=1)

        # Apply random smearing
        smearing = rng.normal(loc=1., scale=electron_smearing_nom_flat)
        

        events = set_ak_column_f32(events, "Electron.pt", ak.unflatten(electron_pt * smearing, arr_shape))

    return events

@electron_smearing_scaling.requires
def electron_smearing_scaling_requires(self: Calibrator, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)
    
@electron_smearing_scaling.setup
def electron_smearing_scaling_setup(
    self: Calibrator,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.electron_scaling_smearing.load(formatter="gzip").decode("utf-8"),
    )
    self.electron_scaling_corrector = correction_set[self.config_inst.x.electron_sf.scale.corrector]
    self.electron_smearing_corrector = correction_set[self.config_inst.x.electron_sf.smearing.corrector]
