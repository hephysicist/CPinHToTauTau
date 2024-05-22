import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={
        "Pileup.nTrueInt"
    },
    produces={
        "pu_weight"
    },
    mc_only=True,
)
def pu_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    nTrueInt = events.Pileup.nTrueInt
    pu_weight = ak.where (self.mc_weight(nTrueInt) != 0,
                          self.data_weight(nTrueInt)/self.mc_weight(nTrueInt) * self.mc2data_norm,
                          0)
    
    events = set_ak_column_f32(events, "pu_weight", pu_weight)
    return events

@pu_weight.setup
def pu_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    """
    Loads the pileup weights added through the requirements and saves them in the
    py:attr:`pu_weight` attribute for simpler access in the actual callable.
    """
    from coffea.lookup_tools import extractor
    ext = extractor()
    data_full_fname = self.config_inst.x.external_files.pileup.data
    data_name = data_full_fname.split('/')[-1].split('.')[0]
    mc_full_fname = self.config_inst.x.external_files.pileup.mc
    mc_name = mc_full_fname.split('/')[-1].split('.')[0]
    ext.add_weight_sets([f'{data_name} pileup {data_full_fname}', f'{mc_name} pileup {mc_full_fname}' ])
    ext.finalize()
    
    self.evaluator = ext.make_evaluator()
    
    mc_integral = 0.
    data_integral = 0.
    for npu in range(0,1000):
        mc_integral += self.evaluator[mc_name](npu)
        data_integral += self.evaluator[data_name](npu)
    
    self.mc_weight = self.evaluator[mc_name]
    self.data_weight = self.evaluator[data_name] 
    self.mc2data_norm = safe_div(mc_integral,data_integral)


@producer(
    uses={"genWeight", optional("LHEWeight.originalXWGTUP")},
    produces={"mc_weight"},
    # only run on mc
    mc_only=True,
)
def get_mc_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Reads the genWeight and LHEWeight columns and makes a decision about which one to save. This
    should have been configured centrally [1] and stored in genWeight, but there are some samples
    where this failed.

    Strategy:

      1. Use LHEWeight.originalXWGTUP when it exists and genWeight is always 1.
      2. In all other cases, use genWeight.

    [1] https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD?rev=99#Weigths
    """
    # determine the mc_weight
    mc_weight = np.sign(events.genWeight)
    if has_ak_column(events, "LHEWeight.originalXWGTUP") and ak.all(events.genWeight == 1.0):
        mc_weight = np.sign(events.LHEWeight.originalXWGTUP)

    # store the column
    events = set_ak_column(events, "mc_weight", mc_weight, value_type=np.float32)

    return events

@producer(
    uses={
        "Muon.pt", "Muon.eta"
    },
    produces={
        "muon_weight"
    },
    mc_only=True,
)
def muon_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    muon_weight = self.muon_sf(self, events.Muon.pt, events.Muon.eta)
    events = set_ak_column(events, "muon_weight", muon_weight, value_type=np.float32)

    return events

@muon_weight.setup
def muon_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    from coffea.lookup_tools import extractor
    ext = extractor()
    full_fname = self.config_inst.x.external_files.muon_correction
    ext.add_weight_sets([f'sf_trig ScaleFactor_trg {full_fname}',
                         f'sf_id ScaleFactor_id {full_fname}',
                         f'sf_iso ScaleFactor_iso {full_fname}'])
    ext.finalize()
    self.evaluator = ext.make_evaluator()
    self.muon_sf = lambda self, pt, eta: self.evaluator['sf_trig'](pt,eta) * self.evaluator['sf_id'](pt,eta) * self.evaluator['sf_iso'](pt,eta)

@producer(
    uses={
        f"Tau.{var}" for var in [
            "pt","eta","decayMode", "genPartFlav"
        ] 
    },
    produces={
        "tau_id_sf"
    },
    mc_only=True,
)
def tau_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Producer for tau scale factors derived by the TAU POG. Requires an external file in the
    config under ``tau_correction``:

        cfg.x.external_files = DotDict.wrap({
            "tau_correction": "/afs/cern.ch/user/s/stzakhar/work/higgs_cp/corrections/tau/POG/TAU/2022_preEE/tau_DeepTau2018v2p5_2022_preEE.json.gz", 
        })

    *get_tau_file* can be adapted in a subclass in case it is stored differently in the external
    files. A correction set named ``"tau_trigger"`` is extracted from it.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun2?rev=113
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/849c6a6efef907f4033715d52290d1a661b7e8f9/POG/TAU
    """
    # the correction tool only supports flat arrays, so convert inputs to flat np view first
    pt = flat_np_view(events.Tau.pt, axis=1)
    abseta = flat_np_view(abs(events.Tau.eta), axis=1)
    dm = flat_np_view(events.Tau.decayMode, axis=1)
    match = flat_np_view(events.Tau.genPartFlav, axis=1)
    
    syst = "nom" # TODO define this systematics inside config file
    deep_tau = self.config_inst.x.deep_tau
    #Get working points of the DeepTau tagger
    #Create scale factor for tau vs jet classifier 
    sf_nom = np.ones_like(pt, dtype=np.float32)
    mask = np.ones_like(pt, dtype=np.float32) #TODO Propagate here the mask for each channel i.e. etau, mutau, tautau 
   
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    #Calculate tau ID scale factors for genuine taus
    # pt, dm, genmatch, jet wp, e wp, syst, sf type
    
       
    tau_mask = flat_np_view(events.Tau.genPartFlav == tau_part_flav["tau_had"])
    sf_nom[tau_mask] = self.id_vs_jet_corrector.evaluate(pt[tau_mask],
                                                         dm[tau_mask],
                                                         match[tau_mask],
                                                         deep_tau.vs_jet,
                                                         deep_tau.vs_e, 
                                                         syst,
                                                         "dm")
    
    #Calculate tau ID scale factors for electron fakes
    #from IPython import embed; embed()
    e_mask = ((match == tau_part_flav["prompt_e"]) | (match == tau_part_flav["tau->e"]))
    
    sf_nom[e_mask] = self.id_vs_e_corrector.evaluate(abseta[e_mask],
                                                     dm[e_mask],
                                                     match[e_mask],
                                                     deep_tau.vs_e, 
                                                     syst)
    
    #Calculate tau ID scale factors for muon fakes
    #from IPython import embed; embed()
    mu_mask = ((match == tau_part_flav["prompt_mu"]) | (match == tau_part_flav["tau->mu"]))
    sf_nom[mu_mask] = self.id_vs_mu_corrector.evaluate(abseta[mu_mask],
                                                      match[mu_mask],
                                                      deep_tau.vs_mu, 
                                                      syst)
    
    #get the original shape of Tau array
    tau_arr_shape = ak.num(events.Tau.pt, axis=1)
    #Unify the shapes of the scale factor array and Tau array 
    sf_shaped = ak.unflatten(sf_nom, tau_arr_shape)
    '''
    If event has more than 1 tau
    '''
    sf_flat = ak.prod(sf_shaped, axis=1)
    #from IPython import embed; embed()
    
    events = set_ak_column(events, "tau_id_sf", sf_flat, value_type=np.float32)
    
    
    return events

@tau_weight.requires
def tau_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@tau_weight.setup
def tau_weight_setup(
    self: Producer,
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
    self.id_vs_jet_corrector    = correction_set[f"{tagger_name}VSjet"]
    self.id_vs_e_corrector      = correction_set[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector     = correction_set[f"{tagger_name}VSmu"]
    self.tes_corrector          = correction_set["tau_energy_scale"]


