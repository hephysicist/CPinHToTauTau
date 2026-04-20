import functools
from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column, has_ak_column, flat_np_view, optional_column as optional
from columnflow.util import maybe_import, load_correction_set

import law
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
    mask = (events.GenZ.pt < 600) | (events.GenZ.mass < 1000) | (events.GenZ.pt > 0) 

    
    zpt = flat_np_view(np.clip(events.GenZ.pt, 0, 600))
    sf_nom = np.ones_like(zpt,dtype=np.float32)
    dataset = self.dataset_inst.name
    if 'DY' in dataset:
        order = self.config_inst.x.dy_ptll_corrs.datasets[dataset]
        sf_nom[mask] = self.zpt_corrector.evaluate(order, zpt[mask], 'nom')
    events = set_ak_column(events, "zpt_weight", sf_nom, value_type=np.float32)
    return events
@zpt_weight.requires
def zpt_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@zpt_weight.setup
def zpt_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    correction_set = load_correction_set(bundle.files.zpt_weight)
    self.zpt_corrector = correction_set[f"DY_pTll_reweighting"]

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
def muon_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@muon_weight.setup
def muon_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
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
def electron_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@electron_weight.setup
def electron_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
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
                    dm_pnet = flat_np_view(tau.decayModePNet, axis=1)
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
                                                       dm_pnet[mask],
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
                    masked_dm_pnet = (dm_pnet == 0) | (dm_pnet == 1) | (dm_pnet == 2) | (dm_pnet == 10) | (dm_pnet == 11)
                    masked_dm = (dm == 0) | (dm == 1) | (dm == 10) | (dm == 11)
                    e_mask = ((genmatch == tau_part_flav["prompt_e"]) | (genmatch == tau_part_flav["tau->e"])) & masked_dm
                    per_ch_sf[e_mask] *= self.id_vs_e_corrector.evaluate(*args_vs_e(e_mask,shift))
                    #Calculate scale factors for tau vs muon classifier 
                    mu_mask = ((genmatch == tau_part_flav["prompt_mu"]) | (genmatch == tau_part_flav["tau->mu"])) & masked_dm
                    per_ch_sf[mu_mask] *= self.id_vs_mu_corrector.evaluate(*args_vs_mu(mu_mask,shift)) 
                    #Calculate tau ID scale factors
                    tau_mask = (genmatch == tau_part_flav["tau_had"]) & masked_dm_pnet
                    per_ch_sf[tau_mask] *= self.id_vs_jet_corrector.evaluate(*args_vs_jet(tau_mask,shift))
                    ch_mask = ak.num(tau, axis=1) > 0
                    shaped_sf = ak.unflatten(per_ch_sf, ak.num(tau.pt, axis=1))
                    sf_values = sf_values * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)       
        events = set_ak_column(events,f"tau_weight_{shift}",sf_values,value_type=np.float32)
                                    
    return events

@tau_weight.requires
def tau_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@tau_weight.setup
def tau_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    cs1 = correctionlib.CorrectionSet.from_string(
        bundle.files.tau_correction.load(formatter="gzip").decode("utf-8"),
    )
    cs2 = correctionlib.CorrectionSet.from_string(
        bundle.files.tau_sf.load(formatter="gzip").decode("utf-8"),
    )
    tagger_name = self.config_inst.x.deep_tau.tagger
    year = self.config_inst.x.year
    tag = self.config_inst.x.long_tag
    self.id_vs_e_corrector      = cs1[f"{tagger_name}VSe"]
    self.id_vs_mu_corrector     = cs1[f"{tagger_name}VSmu"]
    self.id_vs_jet_corrector     = cs2[f"tau_sf_pt-dm_{tagger_name}VSjet_{year}_{tag}"]

@producer(
    uses={
        'event', optional("TauSpinner*") 
    },
    produces={
        "tauspinner_weight"
    },
    mc_only=True,
)
def tauspinner_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    A simple function that sets tauspinner_weight according to the cp_hypothesis
    
    """
    is_signal = ("_htt_" in self.dataset_inst.name)
    if is_signal: 
        proc = self.dataset_inst.processes.get_first().name
        if 'htt_cpo' in proc:
            the_weight = events.TauSpinner.weight_cp_0p5 
        elif 'htt_sm' in proc:
            the_weight = events.TauSpinner.weight_cp_0
        elif 'htt_mm' in proc:
            the_weight = events.TauSpinner.weight_cp_0p25
        elif 'htt_flat' in proc:
            the_weight = np.ones_like(events.event, dtype=np.float32)
        elif (proc == "h_ggf_htt") or (proc == "h_vbf_htt"):
            the_weight = np.ones_like(events.event, dtype=np.float32)
        else: 
            raise NotImplementedError(f'TauSpinner weight for {proc} is not implemented!')
            the_weight = np.ones_like(events.event, dtype=np.float32)
        buf = ak.to_numpy(the_weight)
        if any(np.isnan(buf)):
            warn.warn("tauspinner_weight contains NaNs. Imputing them with zeros.")
            buf[np.isnan(buf)] = 0
            the_weight = buf
    else:
        print("Tauspinner column does not exist for this sample: filling weights with ones")
        the_weight = np.ones_like(events.event, dtype=np.float32)
    events = set_ak_column_f32(events, f"tauspinner_weight", the_weight)
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
    pt = np.clip(flat_np_view(tau.pt, axis=1), 20,300)#fake factors are estimated for this pt region
    dm = flat_np_view(tau.decayModePNet, axis=1)
    n_jets = np.clip(flat_np_view(events.n_jets).copy(),0,2).astype(int) 
    mask = (dm ==0) | (dm==1) | (dm==2) | (dm==10) | (dm==11)
     
    shift_dict = self.config_inst.x.fake_factor_method.shifts
    shifts = ak.cartesian([shift_dict['tau_dm_pnet'],
                           shift_dict['n_j'],
                           shift_dict['shift_name']], axis=0)
   
    args = lambda mask,shift : (pt[mask],
                                dm[mask],
                                n_jets[mask],
                                shift)
    ff_dict = {}
    for the_name in self.config_inst.x.fake_factor_method.columns:
        ff_evaluator = self.fake_factor_qcd if 'qcd' in the_name else self.fake_factor_wjets
        ff_nom = np.zeros_like(pt, dtype=np.float32)
        ff_nom[mask] = ff_evaluator(*args(mask,'nom')) #use these values as nominal
        
        ff_dict[the_name] = ff_nom
        #Evaluate systematics:
        for (the_dm,the_nj,the_syst) in shifts.tolist():
            
            syst_name = '_'.join((the_name,
                            'dm',str(the_dm),
                            'nj',str(the_nj),
                            the_syst))
            
            the_mask = mask & (dm == the_dm) & (n_jets == the_nj)
            ff_syst = ff_nom.copy()
            ff_syst[the_mask] = ff_evaluator(*args(the_mask,the_syst))
            ff_dict[syst_name] = ff_syst
    events = set_ak_column_f32(events, 'ff_weight', ak.zip(ff_dict)) 
    return events

@fake_factors.requires
def fake_factors_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)
    
@fake_factors.setup
def fake_factors_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
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
    
    channel = self.config_inst.channels.names()[0]
    muon = events[f'hcand_{channel}'].lep0 # fake factor method works for nutau/etau channel
    tau = events[f'hcand_{channel}'].lep1 # fake factor method works for nutau/etau channel
    tau_vs_jet_wp = self.config_inst.x.deep_tau.vs_jet[channel]
    
    def prepare(var, low, up, absolute=False):
        if absolute: return np.abs(np.clip(flat_np_view(var), low, up))
        else: return np.clip(flat_np_view(var), low, up)
    
    shifts = ["nom", "up", "down"]
    sf_dict = {}
    for the_shift in shifts:
        mu_syst_name = 'nominal' if the_shift == 'nom' else 'syst' + the_shift
        strig_mu_sf = self.strig_mu.evaluate(prepare(muon.eta,low=-2.4,up=2.4,absolute=True),
                                             prepare(muon.pt,low=26,up=np.inf),
                                             mu_syst_name)
        xtrig_mu_sf = self.xtrig_mu.evaluate(prepare(muon.eta,low=-2.1,up=2.1,absolute=True),
                                             prepare(muon.pt,low=20,up=np.inf),
                                             mu_syst_name)
        xtrig_tau_sf = self.xtrig_tau.evaluate(prepare(tau.pt,low=25,up=np.inf),
                                               flat_np_view(tau.decayMode), #TODO Replace it to PNet for new scale factors: https://gitlab.cern.ch/dwinterb/HiggsDNA/-/tree/master/higgs_dna/systematics/ditau/JSONs/Tau_Trigger/CP
                                               channel,
                                               tau_vs_jet_wp,
                                               'sf',
                                               the_shift)
        sf = np.where(
            flat_np_view(muon.pt) >=26,
            strig_mu_sf,
            xtrig_mu_sf*xtrig_tau_sf)
        events = set_ak_column_f32(events, f"trigger_weight_mutau_{the_shift}",sf) #TODO remve abs when replacuing bugged trigger file
    
    
    return events

@trigger_weight_mutau.requires
def trigger_weight_mutau_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    # Declare external files dependency if not already present
    if "external_files" in reqs:
        return
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@trigger_weight_mutau.setup
def trigger_weight_mutau_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
) -> None:
    # Load and initialize correctionlib objects from external JSON/GZIP files
    bundle = reqs["external_files"]
    import correctionlib
    # Monkey-patch evaluate method to __call__ for convenience
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    muon_file = correctionlib.CorrectionSet.from_file(bundle.files.muon_correction.abspath)
    self.strig_mu = muon_file[self.config_inst.x.muon_sf.trig.corrector]
    
    mutau_file = correctionlib.CorrectionSet.from_file(bundle.files.cross_mutau_mu_leg.abspath)
    self.xtrig_mu = mutau_file[self.config_inst.x.muon_sf.xtrig.corrector]
    
    tau_file  = correctionlib.CorrectionSet.from_file(bundle.files.tau_correction.abspath)
    self.xtrig_tau = tau_file['tau_trigger']

@producer(
    uses={
        'event','hcand_*','n_jets',
    },
    produces={
        'filter_weight'
    },
)
def filter_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    This function applies filter weights to the datasets 
    """
    datasets  = self.dataset_inst.keys
    is_filtered = ['Filtered' in the_name for the_name in datasets]
    if np.any(is_filtered):
        dataset_name = datasets[0].replace('/','')
        the_weight = self.lookup_table[dataset_name]
        print(f'Filter efficiency for {dataset_name} is {the_weight}')
        filter_weight = np.full_like(events.event, the_weight, dtype=np.float32)
    else:
        filter_weight = np.ones_like(events.event, dtype=np.float32)
    events = set_ak_column_f32(events,'filter_weight', filter_weight)
    return events

@filter_weight.requires
def filter_weight_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)

@filter_weight.setup
def filter_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict,
    inputs: dict,
    reader_targets: law.util.InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    self.lookup_table  = bundle.files.filter_eff.load(formatter='yaml')