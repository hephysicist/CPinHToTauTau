import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, DotDict
from law.util import InsertableDict
from columnflow.columnar_util import sorted_indices_from_mask, set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.types import Any
import law

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={f"Jet.{var}" for var in
          ["pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"
           ]},
    produces={
         f"btag_weight_{shift}"
        for shift in ["nom"]
    },
    mc_only=True,
)
def btag_weight(
    self: Producer, 
    events: ak.Array, 
    do_syst: bool,
    **kwargs,  
) -> ak.Array:
    
    shifts = ["central"]
    if do_syst: shifts=[*shifts] 
    sf_values = {}
    tag = self.config_inst.x.tag
    year = self.config_inst.x.year
    btag_wp = self.config_inst.x.btag_working_points[year][tag].deepjet.medium

    #Removing NaNs from discriminat
    dis = events.Jet.btagDeepFlavB 
    nan_mask = np.isnan(dis)
    mask = ~np.isnan(dis)
        
    Jet = events.Jet[mask]

    # base (b-jet) selection
    jet_selections = {
        "jet_pt_20": Jet.pt > 20.0,
        "jet_eta_2.5": abs(Jet.eta) < 2.5,
        "jet_id": Jet.pass_tightID_lep_veto,
        "btag_wp_medium": Jet.btagDeepFlavB >= btag_wp,
    }
    
    jet_obj_mask = ak.ones_like(Jet.pt, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel
    # jet_mask = ((events.Jet.pt >= 20) & 
    #             (abs(events.Jet.eta) < 2.5) & 
    #             (events.Jet.jetId & 0b10 == 0b10))
    

    for the_shift in shifts: sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    # Create sf array template to make copies and dict for finnal results of all systematics
    Jet = Jet[jet_obj_mask]
    flavor = Jet.hadronFlavour
    eta = abs(Jet.eta)
    pt = Jet.pt
    dis = Jet.btagDeepFlavB
    
    for the_shift in shifts: sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)

    #Prepare a tuple with the inputs of the correction evaluator
    btag_sf_args = lambda syst : (syst,flavor,eta,pt,dis)

    #Loop over the shifts and calculate for each shift btag scale factor
    

    for the_shift in shifts:
        sf = ak.ones_like(pt)
        sf = sf * self.btag_sf_corr.evaluate(*btag_sf_args(the_shift))
        sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(sf, 1.)
                    
    rename_systs = {"central" : "nom",} 
             
    for the_shift in shifts: 
        btag_w = sf_values[the_shift]
        w_event = ak.prod(btag_w,axis=1)

        events = set_ak_column_f32(events, f"btag_weight_{rename_systs[the_shift]}", w_event)

    return events

@btag_weight.requires
def btag_weight_requires(
    self: Producer,
    task: law.Task,
    reqs: dict,
    **kwargs,
    ) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

    
@btag_weight.setup
def btag_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
    **kwargs,
) -> None:
   
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.btag_sf_corr.load(formatter="gzip").decode("utf-8"),
    )
    self.btag_sf_corr = correction_set[self.config_inst.x.btag_sf[0]]
