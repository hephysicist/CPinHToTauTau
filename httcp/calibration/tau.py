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

@calibrator(
    uses={f"Tau.{var}" for var in [
                "pt","eta","mass", "decayMode", "genPartFlav"
                ] 
    },
    produces={
        "Tau.pt", "Tau.mass", "Tau.pt_no_tes", "Tau.mass_no_tes" 
    },
    mc_only=True,
)
def tau_energy_scale(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    # fail when running on data
    if self.dataset_inst.is_data:
        raise ValueError("attempt to apply tau energy corrections in data")
    
    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    mass = flat_np_view(events.Tau.mass, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayMode, axis=1)
    match = flat_np_view(events.Tau.genPartFlav, axis=1)
    
    syst = "nom" # TODO define this systematics inside config file
     #Get working points of the DeepTau tagger
    deep_tau = self.config_inst.x.deep_tau
    
    #Create get energy scale correction for each tau
    tes_nom = np.ones_like(pt, dtype=np.float32)
    #Calculate tau ID scale factors for genuine taus
    # pt, eta, dm, genmatch, deep_tau_id, jet_wp, e_wp, syst
    mask2prong = ((dm != 5) & (dm != 6))
    if self.config_inst.x.year in [2022,2023]:
        tes_args = lambda events, mask, deep_tau_obj, syst: (pt[mask],
                                                            abseta[mask],
                                                            dm[mask],
                                                            match[mask],
                                                            deep_tau_obj.tagger,
                                                            deep_tau_obj.vs_jet,
                                                            deep_tau_obj.vs_e,
                                                            syst)
    else:
        tes_args = lambda events, mask, deep_tau_obj, syst: (pt[mask],
                                                            abseta[mask],
                                                            dm[mask],
                                                            match[mask],
                                                            deep_tau_obj.tagger,
                                                            syst)
        
    tes_nom[mask2prong] = self.tes_corrector.evaluate(*tes_args(events, mask2prong,deep_tau, syst))
    
    tes_nom     = np.asarray(tes_nom)
    tau_pt      = np.asarray(ak.flatten(events.Tau.pt))
    tau_mass    = np.asarray(ak.flatten(events.Tau.mass))
    arr_shape = ak.num(events.Tau.pt, axis=1)
    
    events = set_ak_column_f32(events, "Tau.pt_no_tes", ak.unflatten(tau_pt, arr_shape))
    events = set_ak_column_f32(events, "Tau.mass_no_tes", ak.unflatten(tau_mass, arr_shape))
    events = set_ak_column_f32(events, "Tau.pt", ak.unflatten(tau_pt * tes_nom, arr_shape))
    events = set_ak_column_f32(events, "Tau.mass", ak.unflatten(tau_mass * tes_nom, arr_shape))
    return events

@tau_energy_scale.requires
def tau_energy_scale_requires(self: Calibrator, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)
    
@tau_energy_scale.setup
def tau_energy_scale_setup(
    self: Calibrator,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.tau_correction.load(formatter="gzip").decode("utf-8"),
    )
    tagger_name = self.config_inst.x.deep_tau.tagger
    self.tes_corrector = correction_set["tau_energy_scale"]