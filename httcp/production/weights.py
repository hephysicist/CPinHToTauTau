import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask
from httcp.util import compute_eff
import json

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

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

###########################
#### Z pt reweighting #####
###########################
@producer(
    uses={
        "GenZ.*",
    },
    produces={
        "zpt_weight"
    },
    mc_only=True,
)
def zpt_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # is within range
    is_outside_range = (
        ((events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0))
        | ((events.GenZ.pt >= 600.0) | (events.GenZ.mass >= 1000.0))
    )

    sf_nom = ak.ones_like(events.event,dtype=np.float32)
    
    # for safety
    zm  = ak.where(events.GenZ.mass > 1000.0, 999.99, events.GenZ.mass)
    zpt = ak.where(events.GenZ.pt > 600.0, 599.99, events.GenZ.pt)

    processes = self.dataset_inst.processes.names()
    if ak.any(['dy' in proc for proc in processes]):
        sf_nom = ak.where(is_outside_range,
                          1.0,
                          self.zpt_corrector(self,zm,zpt))

    events = set_ak_column(events, "zpt_weight", sf_nom, value_type=np.float32)
    return events

@zpt_weight.setup
def zpt_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    from coffea.lookup_tools import extractor
    ext = extractor()
    full_fname = self.config_inst.x.external_files.zpt_weight
    ext.add_weight_sets([f'zpt_weight zptmass_histo {full_fname}'])
    ext.finalize()
    self.evaluator = ext.make_evaluator()
    self.zpt_corrector = lambda self, mass, pt: self.evaluator['zpt_weight'](mass,pt)
                        

### MUON WEIGHT CALCULATOR ###

@producer(
    uses={
        'hcand_*', 'event'
    },
    produces={
         f"muon_weight_{shift}"
        for shift in ["nom", "up", "down"]
    },
    mc_only=True,
)
def muon_weight(self: Producer, events: ak.Array, do_syst: bool,  **kwargs) -> ak.Array:
    
    shifts = ["nominal"]
    if do_syst: shifts=[*shifts,"systup", "systdown"] 
    sf_values = {}
    for the_shift in shifts: sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        for lep in [field for field in hcand.fields if 'lep' in field]:
            if ch_objects[ch_str][lep] == 'Muon':
                muon = hcand[lep]
                # Create sf array template to make copies and dict for finnal results of all systematics
                pt =  flat_np_view(muon.pt,axis=1) #take the first particle from the hcand pair
                eta =  flat_np_view(abs(muon.eta),axis=1)
                #Prepare a tuple with the inputs of the correction evaluator
                mu_sf_args = lambda syst : (eta,pt,syst)
                #Loop over the shifts and calculate for each shift muon scale factor
                for the_shift in shifts:
                    flat_sf = ak.ones_like(pt)
                    for the_sf in [self.muon_id, self.muon_iso]: #, self.muon_trig]: 
                        flat_sf = flat_sf * the_sf.evaluate(*mu_sf_args(the_shift))
                    shaped_sf = ak.unflatten(flat_sf, ak.num(muon.pt, axis=1))
                    sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)

    rename_systs = {"nominal" : "nom",
                    "systup"  : "up",
                    "systdown": "down"
    }
    for the_shift in shifts: events = set_ak_column_f32(events, f"muon_weight_{rename_systs[the_shift]}", sf_values[the_shift])
    return events

@muon_weight.requires
def muon_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@muon_weight.setup
def muon_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.muon_correction.load(formatter="gzip").decode("utf-8"),
    )
   
    self.muon_id = correction_set[self.config_inst.x.muon_sf.ID.corrector]
    self.muon_iso = correction_set[self.config_inst.x.muon_sf.iso.corrector]
    #self.muon_trig = correction_set[self.config_inst.x.muon_sf.trig.corrector]

# ### ELECTRON WEIGHT CALCULATOR ##
@producer(
    uses={
        'hcand_*', 'event'
    },
    produces={
         f"electron_weight_{shift}"
        for shift in ["nom", "up", "down"]
    },
    mc_only=True,
)
def electron_weight(self: Producer, events: ak.Array, do_syst: bool,  **kwargs) -> ak.Array:
    
    shifts = ["sf"]
    if do_syst: shifts=[*shifts,"sfup", "sfdown"] 
    sf_values = {}
    for the_shift in shifts: sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    year_id = self.config_inst.x.electron_sf.ID.year
    wp_id = self.config_inst.x.electron_sf.ID.wp
    year_trigger = self.config_inst.x.electron_sf.trig.year
    wp_trigger = self.config_inst.x.electron_sf.trig.wp
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        for lep in [field for field in hcand.fields if 'lep' in field]:
            if ch_objects[ch_str][lep] == 'Electron':
                electron = hcand[lep]
                # Create sf array template to make copies and dict for finnal results of all systematics
                pt =  flat_np_view(electron.pt,axis=1) #take the first particle from the hcand pair
                eta =  flat_np_view(electron.eta,axis=1)
                phi =  flat_np_view(electron.phi,axis=1)
                #Prepare a tuple with the inputs of the correction evaluator
                if "2023" in year_id:
                    ele_sf_args_idiso = lambda syst : (year_id,syst,wp_id,eta,pt,phi)
                else:
                    ele_sf_args_idiso = lambda syst :(year_id,syst,wp_id,eta,pt)
                ele_sf_args_trigger = lambda syst :(year_trigger,syst,wp_trigger,eta,pt)
                #Loop over the shifts and calculate for each shift electron scale factor
                for the_shift in shifts:
                    flat_sf = ak.ones_like(pt)
                    flat_sf = flat_sf * self.electron_idiso.evaluate(*ele_sf_args_idiso(the_shift))
                    flat_sf = flat_sf * self.electron_trig.evaluate(*ele_sf_args_trigger(the_shift))
                    shaped_sf = ak.unflatten(flat_sf, ak.num(electron.pt, axis=1))
                    sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)

    rename_systs = {"sf" : "nom",
                    "sfup"  : "up",
                    "sfdown": "down"
    }
    for the_shift in shifts: events = set_ak_column_f32(events, f"electron_weight_{rename_systs[the_shift]}", sf_values[the_shift])
    return events

@electron_weight.requires
def electron_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@electron_weight.setup
def electron_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set_idiso = correctionlib.CorrectionSet.from_string(
        bundle.files.electron_idiso.load(formatter="gzip").decode("utf-8"),
    )
    correction_set_trig = correctionlib.CorrectionSet.from_string(
        bundle.files.electron_trigger.load(formatter="gzip").decode("utf-8"),
    )
   
    self.electron_idiso   = correction_set_idiso[self.config_inst.x.electron_sf.ID.corrector]
    self.electron_trig = correction_set_trig[self.config_inst.x.electron_sf.trig.corrector]


### TAU WEIGHT CALCULATOR ###

@producer(
    uses={
        'event', 'hcand_*',
    },
    produces={
        f"tau_weight_{shift}"
        for shift in ["nom", "up", "down"]
    },
    mc_only=True,
)
def tau_weight(self: Producer, events: ak.Array, do_syst: bool, **kwargs) -> ak.Array:
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

    #Helper function to deal with the case when two taus exist at the events. In that case one should multiply sf values to get the sf per event
    shape_sf = lambda sf: ak.prod(ak.unflatten(sf, 
                                            ak.num(events.Tau.pt, axis=1)), 
                                  axis=1, 
                                  mask_identity=False)
    
    #Make masks for each channel
    shifts = ["nom"]
    if  do_syst: shifts=[*shifts,"up", "down"]         
    sf_values = {}

    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    for shift in shifts:
        sf_values = np.ones_like(events.event, dtype=np.float32)
        for ch_str in channels:
            wp_vs_e   = self.config_inst.x.deep_tau.vs_e.tautau
            wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.tautau
            wp_vs_mu  = self.config_inst.x.deep_tau.vs_mu.tautau
            if ch_str=='mutau':
                wp_vs_e   = self.config_inst.x.deep_tau.vs_e.mutau
                wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.mutau
                wp_vs_mu  = self.config_inst.x.deep_tau.vs_mu.mutau
            elif ch_str=='etau':
                wp_vs_e   = self.config_inst.x.deep_tau.vs_e.etau
                wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.etau
                wp_vs_mu  = self.config_inst.x.deep_tau.vs_mu.etau
            hcand = events[f'hcand_{ch_str}']
            for lep in [field for field in hcand.fields if 'lep' in field]:
                if ch_objects[ch_str][lep]  == 'Tau':
                    tau = hcand[lep]
                    #Prepare flat arrays of the inputs to send into the 
                    pt = flat_np_view(tau.pt, axis=1)
                    eta = flat_np_view(abs(tau.eta), axis=1)
                    dm = flat_np_view(tau.decayMode, axis=1)
                    genmatch = flat_np_view(tau.genPartFlav, axis=1)
                    per_ch_sf = np.ones_like(pt, dtype=np.float32)
                    args_vs_e = lambda mask, syst : (eta[mask],
                                                     dm[mask],
                                                     genmatch[mask],
                                                     wp_vs_e,
                                                     syst)   
                    args_vs_mu = lambda mask, syst : (eta[mask],
                                                      genmatch[mask],
                                                      wp_vs_mu,
                                                      syst)                      
                    args_vs_jet = lambda mask, syst : (pt[mask],
                                                       dm[mask],
                                                       genmatch[mask],
                                                       wp_vs_jet,
                                                       wp_vs_e,
                                                       syst,
                                                       "dm")
                    
                    tau_part_flav = {
                        "prompt_e"  : 1,
                        "prompt_mu" : 2,
                        "tau->e"    : 3,
                        "tau->mu"   : 4,
                        "tau_had"   : 5
                    }
                    #Calculate scale factors for tau vs electron classifier 
                    masked_dm = (dm != 5) & (dm != 6)
                    e_mask = ((genmatch == tau_part_flav["prompt_e"]) | (genmatch == tau_part_flav["tau->e"])) & masked_dm
                    per_ch_sf[e_mask] *= self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask,shift))
                    #Calculate scale factors for tau vs muon classifier 
                    mu_mask = ((genmatch == tau_part_flav["prompt_mu"]) | (genmatch == tau_part_flav["tau->mu"])) & masked_dm
                    per_ch_sf[mu_mask] *= self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask,shift)) 
                    #Calculate tau ID scale factors
                    tau_mask = (genmatch == tau_part_flav["tau_had"]) & masked_dm
                    per_ch_sf[tau_mask] *= self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,shift))

                    ch_mask = ak.num(tau, axis=1) > 0
                    shaped_sf = ak.unflatten(per_ch_sf, ak.num(tau.pt, axis=1))
                    sf_values = sf_values * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)       
        events = set_ak_column(events,f"tau_weight_{shift}",sf_values,value_type=np.float32)
                                    
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

@producer(
    uses={
        'event', optional("TauSpinner*") 
    },
    produces={
        "tauspinner_weight_up", "tauspinner_weight", "tauspinner_weight_down"
    },
    mc_only=True,
)
def tauspinner_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    
    """
    names = ["_up", "", "_down"]
    for the_name in names:
        if  "TauSpinner" in list(events.fields):
            if the_name == "_up": the_weight = events.TauSpinner.weight_cp_0
            elif the_name == "_down": the_weight = events.TauSpinner.weight_cp_0p5
            elif the_name == "":  the_weight =  (events.TauSpinner.weight_cp_0p5 + events.TauSpinner.weight_cp_0)/2.
            else:  raise NotImplementedError('CP hypothesis is not known to the tauspinner weight producer!')   
            buf = ak.to_numpy(the_weight)
            if any(np.isnan(buf)):
                warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
                buf[np.isnan(buf)] = 0
                the_weight = buf
        else:
            print("Tauspinner column does not exist for this sample: filling weights with ones")
            the_weight = np.ones_like(events.event, dtype=np.float32)
        events = set_ak_column_f32(events, f"tauspinner_weight{the_name}", the_weight)
    return events


@producer(
    uses={
        'event','hcand_*','n_jets',
    },
    produces={
        'ff_weight*'
    },
)
def fake_factors(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    
    """
    channel = self.config_inst.channels.names()[0]
    tau = events[f'hcand_{channel}'].lep1 # fake factor method works for nutau/etau channel
    pt = flat_np_view(tau.pt, axis=1)
    dm = flat_np_view(tau.decayModePNet, axis=1)
    n_jets = flat_np_view(events.n_jets).copy()
    n_jets = np.where(n_jets>2,
                      2*np.ones_like(n_jets),
                      n_jets)#scale factors are calculated for nj=0,1,>=2 aka 2
    mask = (pt > 20) & (dm >= 0)
    fake_factors = {}
    for the_name in self.config_inst.x.fake_factor_method.columns:
        for the_shift in self.config_inst.x.fake_factor_method.shifts:
            ff_evaluator = self.fake_factor_qcd if 'qcd' in the_name else self.fake_factor_wjets
            args = lambda mask : (pt[mask],
                                  dm[mask],
                                  n_jets[mask].astype(int),
            )
            ff_vals = np.zeros_like(pt, dtype=np.float32)
            ff_vals[mask] = ff_evaluator(*args(mask))
            events = set_ak_column_f32(events,'_'.join((the_name,the_shift)), ff_vals)
    return events

@fake_factors.requires
def fake_factors_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)
    
@fake_factors.setup
def fake_factors_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    fake_factors = correctionlib.CorrectionSet.from_string(
        bundle.files.fake_factors.load(formatter='json'))
    self.fake_factor_qcd = fake_factors["ff_qcd"]
    self.fake_factor_wjets = fake_factors["ff_wjets"]

### Single Mu or Cross_Mutau SFs CALCULATOR
@producer(
    uses={
        'event', 'hcand_*', 'all_triggers_id', 'triggerID*'
    },
    produces={
        f"trigger_weight_mutau_{the_shift}"
        for the_shift in ["nom", "up", "down"]
    },
    mc_only=True,
)
def trigger_weight_mutau(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    # Initialize dictionaries to hold Data/MC efficiencies for tau and muon legs
    Data_eff_values_tau_xtrig = {}
    MC_eff_values_tau_xtrig   = {}
    Data_eff_values_mu_xtrig  = {}
    MC_eff_values_mu_xtrig    = {}
    Data_eff_values_mu_trig   = {}
    MC_eff_values_mu_trig     = {}

    # Prepare default flat efficiency arrays for each systematic shift
    shifts_tau_xtrig = ["nom", "up", "down"]
    for the_shift in shifts_tau_xtrig:
        Data_eff_values_tau_xtrig[the_shift] = np.ones_like(events.event, dtype=np.float32)
        MC_eff_values_tau_xtrig[the_shift]   = np.ones_like(events.event, dtype=np.float32)
    shifts_mu_xtrig = ["nominal", "systup", "systdown"]
    for the_shift in shifts_mu_xtrig:
        Data_eff_values_mu_xtrig[the_shift] = np.ones_like(events.event, dtype=np.float32)
        MC_eff_values_mu_xtrig[the_shift]   = np.ones_like(events.event, dtype=np.float32)
    shifts_mu_trig = ["nominal", "systup", "systdown"]
    for the_shift in shifts_mu_trig:
        Data_eff_values_mu_trig[the_shift] = np.ones_like(events.event, dtype=np.float32)
        MC_eff_values_mu_trig[the_shift]   = np.ones_like(events.event, dtype=np.float32)

    # Retrieve configured channel names and objects
    channels   = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects

    # -----------------------
    # Tau leg of mutau trigger
    # -----------------------
    for the_shift in shifts_tau_xtrig:
        for ch_str in channels:
            if ch_str == 'mutau':
                # Working point for tau-vs-jet discriminant
                wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.mutau
            hcand = events[f'hcand_{ch_str}']
            # Loop over lepton candidates in the collection
            for lep in [field for field in hcand.fields if 'lep' in field]:
                if ch_objects[ch_str][lep] == 'Tau':
                    # Extract tau candidate
                    tau    = hcand[lep]
                    pt_tau = flat_np_view(tau.pt, axis=1)
                    # For events not passing tau trigger ID 13, use default pt cut
                    pt_tau = ak.where(events.all_triggers_id == 13, pt_tau, 39.598)
                    dm     = flat_np_view(tau.decayMode, axis=1)
                    trigtype = "mutau"

                    # Lambdas to collect args for correction evaluation
                    args_tau_trigger_eff_data = lambda mask, syst: (
                        pt_tau[mask], dm[mask], trigtype, wp_vs_jet, "eff_data", syst
                    )
                    args_tau_trigger_eff_mc = lambda mask, syst: (
                        pt_tau[mask], dm[mask], trigtype, wp_vs_jet, "eff_mc", syst
                    )

                    # Mask out unwanted decay modes
                    masked_dm = (dm != 5) & (dm != 6)
                    # Evaluate Data efficiency
                    flat_eff   = ak.ones_like(pt_tau)
                    flat_eff   = flat_eff * self.id_vs_jet_corrector.evaluate(
                        *args_tau_trigger_eff_data(masked_dm, the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(tau.pt, axis=1))
                    Data_eff_values_tau_xtrig[the_shift] = Data_eff_values_tau_xtrig[the_shift]*shaped_eff
                    # Evaluate MC efficiency
                    flat_eff   = ak.ones_like(pt_tau)
                    flat_eff   = flat_eff * self.id_vs_jet_corrector.evaluate(
                        *args_tau_trigger_eff_mc(masked_dm, the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(tau.pt, axis=1))
                    MC_eff_values_tau_xtrig[the_shift] = MC_eff_values_tau_xtrig[the_shift]*shaped_eff

    # ---------------------------------
    # Muon of Single mu trigger 
    # ---------------------------------
    for the_shift in shifts_mu_trig:
        for ch_str in channels:
            if ch_str == 'mutau':
                wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.mutau
            hcand = events[f'hcand_{ch_str}']
            for lep in [field for field in hcand.fields if 'lep' in field]:
                if ch_objects[ch_str][lep] == 'Muon':
                    muon = hcand[lep]
                    pt   = flat_np_view(muon.pt, axis=1)
                    # Default pt cut for events without triggerID 132
                    pt   = ak.where(events.all_triggers_id == 132, pt, 26)
                    pt   = flat_np_view(pt)
                    eta  = flat_np_view(abs(muon.eta), axis=1)
                    # Lambda for muon trigger SF args
                    muon_sf_args_trigger_eff = lambda syst: (eta, pt, syst)

                    # Data efficiency
                    flat_eff   = ak.ones_like(pt)
                    flat_eff   = flat_eff * self.muon_trig_data_eff.evaluate(
                        *muon_sf_args_trigger_eff(the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))
                    Data_eff_values_mu_trig[the_shift] = Data_eff_values_mu_trig[the_shift]*shaped_eff
                    # MC efficiency
                    flat_eff   = ak.ones_like(pt)
                    flat_eff   = flat_eff * self.muon_trig_mc_eff.evaluate(
                        *muon_sf_args_trigger_eff(the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))
                    MC_eff_values_mu_trig[the_shift] = MC_eff_values_mu_trig[the_shift]*shaped_eff

    # ---------------------------------
    # Muon leg of mutau trigger 
    # Cross mu leg SFs
    # ---------------------------------
    for the_shift in shifts_mu_xtrig:
        for ch_str in channels:
            if ch_str == 'mutau':
                wp_vs_jet = self.config_inst.x.deep_tau.vs_jet.mutau
            hcand = events[f'hcand_{ch_str}']
            for lep in [field for field in hcand.fields if 'lep' in field]:
                if ch_objects[ch_str][lep] == 'Muon':
                    muon = hcand[lep]
                    pt   = flat_np_view(muon.pt, axis=1)
                    # Default pt and eta for cross trigger SF
                    pt   = ak.where(events.all_triggers_id == 13, pt, 21)
                    eta  = flat_np_view(abs(muon.eta), axis=1)
                    eta  = ak.where(events.all_triggers_id == 13, eta, 0)
                    muon_sf_args_trigger_eff = lambda syst: (eta, pt, syst)

                    # Data efficiency
                    flat_eff   = ak.ones_like(pt)
                    flat_eff   = flat_eff * self.mutau_muleg_Data_eff.evaluate(
                        *muon_sf_args_trigger_eff(the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))
                    Data_eff_values_mu_xtrig[the_shift] = Data_eff_values_mu_xtrig[the_shift]*shaped_eff
                    # MC efficiency
                    flat_eff   = ak.ones_like(pt)
                    flat_eff   = flat_eff * self.mutau_muleg_MC_eff.evaluate(
                        *muon_sf_args_trigger_eff(the_shift)
                    )
                    shaped_eff = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))
                    MC_eff_values_mu_xtrig[the_shift] = MC_eff_values_mu_xtrig[the_shift]*shaped_eff

    # Mapping of systematic labels
    rename_systs = {"nominal": "nom", "systup": "up", "systdown": "down"}
    # Write out intermediate efficiency columns into events
    for the_shift in shifts_tau_xtrig:
        events = set_ak_column_f32(events, f"tau_xtrig_Data_eff_{the_shift}", Data_eff_values_tau_xtrig[the_shift])
        events = set_ak_column_f32(events, f"tau_xtrig_MC_eff_{the_shift}",   MC_eff_values_tau_xtrig[the_shift])
    for the_shift in shifts_mu_xtrig:
        short = rename_systs[the_shift]
        events = set_ak_column_f32(events, f"muon_xtrig_Data_eff_{short}", Data_eff_values_mu_xtrig[the_shift])
        events = set_ak_column_f32(events, f"muon_xtrig_MC_eff_{short}",   MC_eff_values_mu_xtrig[the_shift])
    for the_shift in shifts_mu_trig:
        short = rename_systs[the_shift]
        events = set_ak_column_f32(events, f"muon_trig_Data_eff_{short}", Data_eff_values_mu_trig[the_shift])
        events = set_ak_column_f32(events, f"muon_trig_MC_eff_{short}",   MC_eff_values_mu_trig[the_shift])

    # Define triggers passed masks
    Pass_mu_trig    = ak.where(events.triggerID_mu > 0   , 1, 0)
    Pass_mutau_trig = ak.where(events.triggerID_mutau > 0, 1, 0)

    # Organize inputs for final SF calculation
    eff_inputs = {
        'data': {
            'nom':  (events.muon_trig_Data_eff_nom,  events.muon_xtrig_Data_eff_nom,  events.tau_xtrig_Data_eff_nom),
            'up':   (events.muon_trig_Data_eff_up,   events.muon_xtrig_Data_eff_up,   events.tau_xtrig_Data_eff_up),
            'down': (events.muon_trig_Data_eff_down, events.muon_xtrig_Data_eff_down, events.tau_xtrig_Data_eff_down),
        },
        'mc': {
            'nom':  (events.muon_trig_MC_eff_nom,    events.muon_xtrig_MC_eff_nom,    events.tau_xtrig_MC_eff_nom),
            'up':   (events.muon_trig_MC_eff_up,     events.muon_xtrig_MC_eff_up,     events.tau_xtrig_MC_eff_up),
            'down': (events.muon_trig_MC_eff_down,   events.muon_xtrig_MC_eff_down,   events.tau_xtrig_MC_eff_down),
        }
    }

    # # Compute final scale factors and add to events
    # for var in ('nom', 'up', 'down'):
    #     eff_data = compute_eff(Pass_mu_trig, Pass_mutau_trig, *eff_inputs['data'][var])
    #     eff_mc   = compute_eff(Pass_mu_trig, Pass_mutau_trig, *eff_inputs['mc'][var])
    #     sf       = eff_data / eff_mc
    #     events   = set_ak_column_f32(events, f"trigger_weight_mutau_{var}", sf)
    
    # NOM
    # data 
    eff_trig_mu_Data_nom = events.muon_trig_Data_eff_nom
    eff_xtrig_mu_Data_nom = events.muon_xtrig_Data_eff_nom
    eff_xtrig_tau_Data_nom = events.tau_xtrig_Data_eff_nom
    # MC
    eff_trig_mu_MC_nom = events.muon_trig_MC_eff_nom
    eff_xtrig_mu_MC_nom = events.muon_xtrig_MC_eff_nom
    eff_xtrig_tau_MC_nom = events.tau_xtrig_MC_eff_nom
    # UP
    # data 
    eff_trig_mu_Data_up = events.muon_trig_Data_eff_up
    eff_xtrig_mu_Data_up = events.muon_xtrig_Data_eff_up
    eff_xtrig_tau_Data_up = events.tau_xtrig_Data_eff_up
    # MC
    eff_trig_mu_MC_up = events.muon_trig_MC_eff_up
    eff_xtrig_mu_MC_up = events.muon_xtrig_MC_eff_up
    eff_xtrig_tau_MC_up = events.tau_xtrig_MC_eff_up
    # DOWN
    # data
    eff_trig_mu_Data_down = events.muon_trig_Data_eff_down
    eff_xtrig_mu_Data_down = events.muon_xtrig_Data_eff_down
    eff_xtrig_tau_Data_down = events.tau_xtrig_Data_eff_down
    # MC
    eff_trig_mu_MC_down = events.muon_trig_MC_eff_down
    eff_xtrig_mu_MC_down = events.muon_xtrig_MC_eff_down
    eff_xtrig_tau_MC_down = events.tau_xtrig_MC_eff_down
    
    # If  eff_trig_mu_Data_nom < eff_xtrig_mu_Data_nom and  eff_trig_mu_MC_nom <  eff_xtrig_mu_MC_nom
    
    one_like_array = np.ones_like(eff_trig_mu_Data_nom, dtype=np.float32)
    
    # Evauation of the nominal efficiencies for data and mc
    eff_data = Pass_mu_trig*eff_trig_mu_Data_nom - Pass_mu_trig*Pass_mutau_trig*eff_trig_mu_Data_nom*eff_xtrig_tau_Data_nom + Pass_mutau_trig*eff_xtrig_mu_Data_nom*eff_xtrig_tau_Data_nom

    eff_mc = Pass_mu_trig*eff_trig_mu_MC_nom - Pass_mu_trig*Pass_mutau_trig*eff_trig_mu_MC_nom*eff_xtrig_tau_MC_nom + Pass_mutau_trig*eff_xtrig_mu_MC_nom*eff_xtrig_tau_MC_nom
    
    SF = eff_data/eff_mc
    # Evauation of the partial derivaties for data and mc
    delta_eff_data_nom_delta_eff_trig_mu = Pass_mu_trig - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_tau_Data_nom
    delta_eff_data_nom_delta_eff_xtrig_mu = Pass_mutau_trig*eff_xtrig_tau_Data_nom
    delta_eff_data_nom_delta_eff_xtrig_tau = Pass_mutau_trig*eff_xtrig_mu_Data_nom - Pass_mu_trig*Pass_mutau_trig*eff_trig_mu_Data_nom 
    
    delta_eff_data_mu = (eff_trig_mu_Data_up - eff_trig_mu_Data_nom)
    delta_eff_data_xmu = (eff_xtrig_mu_Data_up - eff_xtrig_mu_Data_nom)
    delta_eff_data_xtau =  (eff_xtrig_tau_Data_up - eff_xtrig_tau_Data_nom)

    delta_eff_mc_nom_delta_eff_trig_mu = Pass_mu_trig - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_tau_MC_nom
    delta_eff_mc_nom_delta_eff_xtrig_mu = Pass_mutau_trig*eff_xtrig_tau_MC_nom
    delta_eff_mc_nom_delta_eff_xtrig_tau = Pass_mutau_trig*eff_xtrig_mu_MC_nom - Pass_mu_trig*Pass_mutau_trig*eff_trig_mu_MC_nom
    
    delta_eff_mc_mu   = (eff_trig_mu_MC_up - eff_trig_mu_MC_nom)
    delta_eff_mc_xmu  = (eff_xtrig_mu_MC_up - eff_xtrig_mu_MC_nom)
    delta_eff_mc_xtau =  (eff_xtrig_tau_MC_up - eff_xtrig_tau_MC_nom)
    
    # Evauation of the partial derivaties for data and mc
    Delta_eff_data = delta_eff_data_nom_delta_eff_trig_mu*delta_eff_data_mu + delta_eff_data_nom_delta_eff_xtrig_mu*delta_eff_data_xmu + delta_eff_data_nom_delta_eff_xtrig_tau*delta_eff_data_xtau
    Delta_eff_mc   = delta_eff_mc_nom_delta_eff_trig_mu*delta_eff_mc_mu + delta_eff_mc_nom_delta_eff_xtrig_mu*delta_eff_mc_xmu + delta_eff_mc_nom_delta_eff_xtrig_tau*delta_eff_mc_xtau
    
    DELTA_EFF_sq = np.sqrt((Delta_eff_data/eff_data)**2 + (Delta_eff_mc/eff_mc)**2)
    DELTA_EFF    = SF*DELTA_EFF_sq
    
    SF_UP   = SF + DELTA_EFF
    SF_DOWN = SF - DELTA_EFF
    
    # If  eff_mu_trig_Data_nom >  eff_mu_xtrig_Data_nom and  eff_mu_trig_MC_nom > eff_mu_xtrig_MC_nom

    # Evauation of the nominal efficiencies for data and mc

    eff_data_1 = Pass_mu_trig*eff_trig_mu_Data_nom - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_mu_Data_nom*eff_xtrig_tau_Data_nom + Pass_mutau_trig*eff_xtrig_mu_Data_nom*eff_xtrig_tau_Data_nom

    eff_mc_1 = Pass_mu_trig*eff_trig_mu_MC_nom - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_mu_MC_nom*eff_xtrig_tau_MC_nom + Pass_mutau_trig*eff_xtrig_mu_MC_nom*eff_xtrig_tau_MC_nom

    SF_1 = eff_data_1/eff_mc_1
    
    # Evauation of the partial derivaties for data and mc
    delta_eff_data_nom_delta_eff_trig_mu = Pass_mu_trig 
    delta_eff_data_nom_delta_eff_xtrig_mu = Pass_mutau_trig*eff_xtrig_tau_Data_nom - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_tau_Data_nom
    delta_eff_data_nom_delta_eff_xtrig_tau = Pass_mutau_trig*eff_xtrig_mu_Data_nom - Pass_mutau_trig*eff_xtrig_mu_Data_nom

    delta_eff_data_mu = (eff_trig_mu_Data_up - eff_trig_mu_Data_nom)
    delta_eff_data_xmu = (eff_xtrig_mu_Data_up - eff_xtrig_mu_Data_nom)
    delta_eff_data_xtau =  (eff_xtrig_tau_Data_up - eff_xtrig_tau_Data_nom)

    delta_eff_mc_nom_delta_eff_trig_mu = Pass_mu_trig 
    delta_eff_mc_nom_delta_eff_xtrig_mu = Pass_mutau_trig*eff_xtrig_tau_MC_nom - Pass_mu_trig*Pass_mutau_trig*eff_xtrig_tau_MC_nom
    delta_eff_mc_nom_delta_eff_xtrig_tau = Pass_mutau_trig*eff_xtrig_mu_MC_nom - Pass_mutau_trig*eff_xtrig_mu_MC_nom

    delta_eff_mc_mu   = (eff_trig_mu_MC_up - eff_trig_mu_MC_nom)
    delta_eff_mc_xmu  = (eff_xtrig_mu_MC_up - eff_xtrig_mu_MC_nom)
    delta_eff_mc_xtau =  (eff_xtrig_tau_MC_up - eff_xtrig_tau_MC_nom)

    # Evauation of the partial derivaties for data and mc
    Delta_eff_data = delta_eff_data_nom_delta_eff_trig_mu*delta_eff_data_mu + delta_eff_data_nom_delta_eff_xtrig_mu*delta_eff_data_xmu + delta_eff_data_nom_delta_eff_xtrig_tau*delta_eff_data_xtau
    Delta_eff_mc   = delta_eff_mc_nom_delta_eff_trig_mu*delta_eff_mc_mu + delta_eff_mc_nom_delta_eff_xtrig_mu*delta_eff_mc_xmu + delta_eff_mc_nom_delta_eff_xtrig_tau*delta_eff_mc_xtau

    DELTA_EFF_sq = np.sqrt((Delta_eff_data/eff_data_1)**2 + (Delta_eff_mc/eff_mc_1)**2)
    DELTA_EFF    = SF_1*DELTA_EFF_sq
    
    SF_1_UP   = SF_1 + DELTA_EFF
    SF_1_DOWN = SF_1 - DELTA_EFF
    
    mask = (eff_trig_mu_Data_nom < eff_xtrig_mu_Data_nom)
    SF_final = ak.where(mask,SF,SF_1)
    SF_final_UP = ak.where(mask,SF_UP,SF_1_UP)
    SF_final_DOWN = ak.where(mask,SF_DOWN,SF_1_DOWN)

    events = set_ak_column_f32(events, f"trigger_weight_mutau_nom", SF_final)
    events = set_ak_column_f32(events, f"trigger_weight_mutau_up", SF_final_UP)
    events = set_ak_column_f32(events, f"trigger_weight_mutau_down", SF_final_DOWN)
    
    return events

@trigger_weight_mutau.requires
def trigger_weight_mutau_requires(self: Producer, reqs: dict) -> None:
    # Declare external files dependency if not already present
    if "external_files" in reqs:
        return
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@trigger_weight_mutau.setup
def trigger_weight_mutau_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    # Load and initialize correctionlib objects from external JSON/GZIP files
    bundle = reqs["external_files"]
    import correctionlib
    # Monkey-patch evaluate method to __call__ for convenience
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    # 1) Single muon trigger corrections
    muon_corr  = bundle.files.HLT_mu_eff.load(formatter="json")
    muon_json  = json.dumps(muon_corr)
    cs_mu      = correctionlib.CorrectionSet.from_string(muon_json)
    self.muon_trig_data_eff = cs_mu[self.config_inst.x.muon_sf.trig_data_eff.corrector]
    self.muon_trig_mc_eff   = cs_mu[self.config_inst.x.muon_sf.trig_mc_eff.corrector]

    # 2) Muon leg of mutau trigger corrections
    muon_corr2 = bundle.files.cross_mutau_mu_leg.load(formatter="json")
    muon_json2 = json.dumps(muon_corr2)
    cs_mu2     = correctionlib.CorrectionSet.from_string(muon_json2)
    self.mutau_muleg_Data_eff = cs_mu2[self.config_inst.x.muon_sf.Data_eff_mutau.corrector]
    self.mutau_muleg_MC_eff   = cs_mu2[self.config_inst.x.muon_sf.MC_eff_mutau.corrector]

    # 3) Tau leg of mutau trigger corrections
    tau_corr_bytes   = bundle.files.tau_correction.load(formatter="gzip")
    tau_corr_str     = tau_corr_bytes.decode("utf-8")
    cs_tau           = correctionlib.CorrectionSet.from_string(tau_corr_str)
    self.id_vs_jet_corrector = cs_tau["tau_trigger"]


