import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask
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

    [1] 
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD?rev=99#Weigths

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
        'hcand_*', 'event', 'all_triggers_id',
    },
    produces={
        name
        for shift in ["nom", "up", "down"]
        for name in (
            f"muon_weight_{shift}",
            f"muon_Data_eff_{shift}",
            f"muon_MC_eff_{shift}",
        )
    },
    mc_only=True,
)
def muon_weight(self: Producer, events: ak.Array, do_syst: bool,  **kwargs) -> ak.Array:
    
    shifts = ["nominal"]
    if do_syst: shifts=[*shifts,"systup", "systdown"] 
    shifts_eff = ["nominal"]
    if do_syst: shifts_eff=[*shifts_eff,"systup", "systdown"] 
    sf_values = {}
    Data_eff_values = {}
    MC_eff_values = {}
    for the_shift in shifts: 
        sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    for the_shift in shifts_eff: 
        Data_eff_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
        MC_eff_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        for lep in [field for field in hcand.fields if 'lep' in field]:
            if ch_objects[ch_str][lep] == 'Muon':
                muon = hcand[lep]
                # Create sf array template to make copies and dict for finnal results of all systematics
                pt =  flat_np_view(muon.pt,axis=1)
                pt = ak.where(events.all_triggers_id == 132, pt, 26)
                pt = flat_np_view(pt)
                eta =  flat_np_view(abs(muon.eta),axis=1) 
                O_mu = (pt > 26) & (abs(eta)<2.4)
                #Prepare a tuple with the inputs of the correction evaluator
                mu_sf_args = lambda syst : (eta,pt,syst)
                
                muon_sf_args_trigger_eff = lambda syst : (eta,pt,syst)
                #Loop over the shifts and calculate for each shift muon scale factor
                for the_shift in shifts:
                    flat_sf = ak.ones_like(pt)
                    for the_sf in [self.muon_id, self.muon_iso]: # self.muon_trig]: 
                        flat_sf = flat_sf * the_sf.evaluate(*mu_sf_args(the_shift))
                    shaped_sf = ak.unflatten(flat_sf, ak.num(muon.pt, axis=1))
                    sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)
                    sf_values[the_shift] =  ak.where(events.all_triggers_id == 132, sf_values[the_shift],1.)
                for the_shift in shifts_eff:
                    flat_eff              = ak.ones_like(pt)
                    flat_eff              = flat_eff * self.muon_trig_Data_eff.evaluate(*muon_sf_args_trigger_eff(the_shift))
                    shaped_eff            = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))                    
                    Data_eff_values[the_shift] = Data_eff_values[the_shift] * ak.fill_none(ak.firsts(shaped_eff,axis=1), 0)
                    Data_eff_values[the_shift] =  ak.where(O_mu, Data_eff_values[the_shift],0)
                    #evaluating the trigger MC efficiency
                    flat_eff              = ak.ones_like(pt)
                    flat_eff              = flat_eff * self.muon_trig_MC_eff.evaluate(*muon_sf_args_trigger_eff(the_shift))
                    shaped_eff            = ak.unflatten(flat_eff, ak.num(muon.pt, axis=1))                    
                    MC_eff_values[the_shift] = MC_eff_values[the_shift] * ak.fill_none(ak.firsts(shaped_eff,axis=1), 0)
                    MC_eff_values[the_shift] =  ak.where(O_mu, MC_eff_values[the_shift],0)  
    rename_systs = {"nominal" : "nom",
                    "systup"  : "up",
                    "systdown": "down"
    }             
    for the_shift in shifts: 
        events = set_ak_column_f32(events, f"muon_weight_{rename_systs[the_shift]}", sf_values[the_shift])
    for the_shift in shifts_eff: 
        events = set_ak_column_f32(events, f"muon_Data_eff_{rename_systs[the_shift]}", Data_eff_values[the_shift])
        events = set_ak_column_f32(events, f"muon_MC_eff_{rename_systs[the_shift]}", MC_eff_values[the_shift])
    #Setting up and down variations, only stat uncertainties for now
    events = set_ak_column_f32(events, f"muon_Data_eff_up", Data_eff_values["systup"])
    events = set_ak_column_f32(events, f"muon_Data_eff_down", Data_eff_values["systdown"])
    events = set_ak_column_f32(events, f"muon_MC_eff_up", MC_eff_values["systup"])
    events = set_ak_column_f32(events, f"muon_MC_eff_down", MC_eff_values["systdown"])


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

    # 1) load the JSON file into a Python dict
    muon_corr = bundle.files.HLT_mu_eff.load(formatter="json")

    # 2) serialize the dict back to a JSON string
    muon_json = json.dumps(muon_corr)

    # 3) build your CorrectionSet
    correction_set_trig = correctionlib.CorrectionSet.from_string(muon_json)
    self.muon_trig_Data_eff = correction_set_trig[self.config_inst.x.muon_sf.trig_data_eff.corrector]
    self.muon_trig_MC_eff = correction_set_trig[self.config_inst.x.muon_sf.trig_mc_eff.corrector]

# ### ELECTRON WEIGHT CALCULATOR ##
@producer(
    uses={
        'hcand_*', 'event'
    },
    produces={
        name
        for shift in ["nom", "up", "down"]
        for name in (
            f"electron_weight_{shift}",
            f"electron_Data_eff_{shift}",
            f"electron_MC_eff_{shift}",
        )
    },
    mc_only=True,
)
def electron_weight(self: Producer, events: ak.Array, do_syst: bool,  **kwargs) -> ak.Array:
    
    shifts = ["sf"]
    if do_syst: shifts=[*shifts,"sfup", "sfdown"]
    shifts_eff = ["nom"]
    if do_syst: shifts_eff=[*shifts_eff,"up", "down"] 
    sf_values = {}
    Data_eff_values = {}
    MC_eff_values = {}
    for the_shift in shifts: 
        sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    for the_shift in shifts_eff: 
        Data_eff_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
        MC_eff_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    year_id = self.config_inst.x.electron_sf.ID.year
    wp_id = self.config_inst.x.electron_sf.ID.wp
    wp_id_RecoBelow20 = "RecoBelow20"
    wp_id_Reco20to75 = "Reco20to75"
    wp_id_RecoAbove75 = "RecoAbove75"
    year_trigger = self.config_inst.x.electron_sf.trig.year
    wp_trigger = self.config_inst.x.electron_sf.trig.wp
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        for lep in [field for field in hcand.fields if 'lep' in field]:
            if ch_objects[ch_str][lep] == 'Electron':
                electron = hcand[lep]
                # Create sf array template to make copies and dict for finnal results of all systematics
                pt            = flat_np_view(electron.pt,axis=1) #take the first particle from the hcand pair
                eta           = flat_np_view(electron.eta,axis=1)
                phi           = flat_np_view(electron.phi,axis=1)
                pt_20_mask    = flat_np_view(electron.pt < 20)
                pt_20_75_mask = flat_np_view((electron.pt >= 20) & (electron.pt < 75))
                pt_75_mask    = flat_np_view(electron.pt >= 75)
                pt_20         = flat_np_view(ak.where(pt_20_mask,pt,19))
                pt_20_75      = flat_np_view(ak.where(pt_20_75_mask,pt,21))
                pt_75         = flat_np_view(ak.where(pt_75_mask, pt, 76))
                pt_trigg      = ak.where(events.all_triggers_id == 111, pt, 26)
                pt_trigg      = flat_np_view(pt_trigg)
                O_e = (pt > 31) & (abs(eta)<2.4)
                
                #Prepare a tuple with the inputs of the correction evaluator
                if '2023' in year_id :
                    ele_sf_args_idiso       = lambda syst : (year_id,syst,wp_id,eta,pt,phi)
                    ele_sf_args_RecoBelow20 = lambda syst : (year_id,syst,wp_id_RecoBelow20,eta,pt_20,phi)
                    ele_sf_args_Reco20to75  = lambda syst : (year_id,syst,wp_id_Reco20to75,eta,pt_20_75,phi)
                    ele_sf_args_RecoAbove75 = lambda syst : (year_id,syst,wp_id_RecoAbove75,eta,pt_75,phi)   
                else :
                    ele_sf_args_idiso       = lambda syst : (year_id,syst,wp_id,eta,pt)
                    ele_sf_args_RecoBelow20 = lambda syst : (year_id,syst,wp_id_RecoBelow20,eta,pt_20)
                    ele_sf_args_Reco20to75  = lambda syst : (year_id,syst,wp_id_Reco20to75,eta,pt_20_75)
                    ele_sf_args_RecoAbove75 = lambda syst : (year_id,syst,wp_id_RecoAbove75,eta,pt_75)

                ele_sf_args_trigger = lambda syst : (year_trigger,syst,wp_trigger,eta,pt_trigg)
                ele_sf_args_trigger_eff = lambda syst : (year_trigger,syst,wp_trigger,eta,pt_trigg)

                #Loop over the shifts and calculate for each shift electron scale factor
                for the_shift in shifts:
                    #Extracting the SFs
                    flat_sf              = ak.ones_like(pt)
                    flat_sf              = flat_sf * self.electron_idiso.evaluate(*ele_sf_args_idiso(the_shift))
                    RecoBelow20_sf       = self.electron_idiso.evaluate(*ele_sf_args_RecoBelow20(the_shift))
                    Reco20to75_sf        = self.electron_idiso.evaluate(*ele_sf_args_Reco20to75(the_shift))
                    RecoAbove75_sf       = self.electron_idiso.evaluate(*ele_sf_args_RecoAbove75(the_shift))
                    flat_sf              = ak.where(pt_20_mask,flat_sf*RecoBelow20_sf,flat_sf)
                    flat_sf              = ak.where(pt_20_75_mask,flat_sf*Reco20to75_sf,flat_sf)
                    flat_sf              = ak.where(pt_75_mask,flat_sf*RecoAbove75_sf,flat_sf)
                    #flat_sf              = flat_sf * self.electron_trig.evaluate(*ele_sf_args_trigger(the_shift))
                    shaped_sf            = ak.unflatten(flat_sf, ak.num(electron.pt, axis=1))
                    sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(ak.firsts(shaped_sf,axis=1), 1.)
                    sf_values[the_shift] =  ak.where(events.all_triggers_id == 111, sf_values[the_shift],1)
                    #evaluating the trigger Data efficiency
                for the_shift in shifts_eff:
                    flat_eff              = ak.ones_like(pt)
                    flat_eff              = flat_eff * self.electron_trig_Data_eff.evaluate(*ele_sf_args_trigger_eff(the_shift))
                    shaped_eff            = ak.unflatten(flat_eff, ak.num(electron.pt, axis=1))                    
                    Data_eff_values[the_shift] = Data_eff_values[the_shift] * ak.fill_none(ak.firsts(shaped_eff,axis=1), 0)
                    Data_eff_values[the_shift] =  ak.where(O_e, Data_eff_values[the_shift],0)
                    #evaluating the trigger MC efficiency
                    flat_eff              = ak.ones_like(pt)
                    flat_eff              = flat_eff * self.electron_trig_MC_eff.evaluate(*ele_sf_args_trigger_eff(the_shift))
                    shaped_eff            = ak.unflatten(flat_eff, ak.num(electron.pt, axis=1))                    
                    MC_eff_values[the_shift] = MC_eff_values[the_shift] * ak.fill_none(ak.firsts(shaped_eff,axis=1), 0)
                    MC_eff_values[the_shift] =  ak.where(O_e, MC_eff_values[the_shift],0)
    
    rename_systs = {"sf" : "nom",
                    "sfup"  : "up",
                    "sfdown": "down"
    }
    for the_shift in shifts: 
        events = set_ak_column_f32(events, f"electron_weight_{rename_systs[the_shift]}", sf_values[the_shift])
    for the_shift in shifts_eff: 
        events = set_ak_column_f32(events, f"electron_Data_eff_{the_shift}", Data_eff_values[the_shift])
    for the_shift in shifts_eff: 
        events = set_ak_column_f32(events, f"electron_MC_eff_{the_shift}", MC_eff_values[the_shift])
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
    #self.electron_trig = correction_set_trig[self.config_inst.x.electron_sf.trig.corrector]
    self.electron_trig_Data_eff = correction_set_trig[self.config_inst.x.electron_sf.Data_eff.corrector]
    self.electron_trig_MC_eff = correction_set_trig[self.config_inst.x.electron_sf.MC_eff.corrector]

################################
#### Trigger SF evaluator ######
################################
@producer(
    uses={
        'event', 'triggerID*'},
    produces={
         f"Trigger_SF_{shift}"
        for shift in ["nom", "up", "down"]
    },
    mc_only=True,
)
def trigger_sf(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    # Nominal
    muon_Data_eff_nom  = events.muon_Data_eff_nom
    electron_Data_eff_nom = events.electron_Data_eff_nom
    muon_MC_eff_nom  = events.muon_MC_eff_nom
    electron_MC_eff_nom = events.electron_MC_eff_nom
    # Up
    muon_Data_eff_up  = events.muon_Data_eff_up
    electron_Data_eff_up = events.electron_Data_eff_up
    muon_MC_eff_up  = events.muon_MC_eff_up
    electron_MC_eff_up = events.electron_MC_eff_up
    # Down
    muon_Data_eff_down  = events.muon_Data_eff_down
    electron_Data_eff_down = events.electron_Data_eff_down
    muon_MC_eff_down  = events.muon_MC_eff_down
    electron_MC_eff_down = events.electron_MC_eff_down
    
    # Define triggers passed masks
    Pass_mu_trig = ak.where(events.triggerID_mu > 0, 1, 0)
    Pass_e_trig  = ak.where(events.triggerID_e >  0, 1, 0)
    
    # SF evaluation
    # nom 
    eff_data  = Pass_mu_trig*muon_Data_eff_nom + Pass_e_trig*electron_Data_eff_nom - Pass_mu_trig*muon_Data_eff_nom*Pass_e_trig*electron_Data_eff_nom
    eff_MC  = Pass_mu_trig*muon_MC_eff_nom + Pass_e_trig*electron_MC_eff_nom - Pass_mu_trig*muon_MC_eff_nom*Pass_e_trig*electron_MC_eff_nom
    
    SF_nom =  eff_data/eff_MC
    
    deff_data_delectron_Data_eff_nom = Pass_e_trig  - Pass_e_trig*Pass_mu_trig*muon_Data_eff_nom
    deff_data_dmuon_Data_eff_nom     = Pass_mu_trig - Pass_e_trig*Pass_mu_trig*electron_Data_eff_nom
    
    delta_el_data = (electron_Data_eff_up - electron_Data_eff_nom)
    delta_mu_data = (muon_Data_eff_up - muon_Data_eff_nom)
    
    deff_data = deff_data_delectron_Data_eff_nom*delta_el_data + deff_data_dmuon_Data_eff_nom*delta_mu_data
    
    deff_data_rel = deff_data/eff_data 
    
    deff_MC_delectron_MC_eff_nom = Pass_e_trig  - Pass_e_trig*Pass_mu_trig*muon_MC_eff_nom
    deff_MC_dmuon_MC_eff_nom     = Pass_mu_trig - Pass_e_trig*Pass_mu_trig*electron_MC_eff_nom
    
    delta_el_MC = (electron_MC_eff_up - electron_MC_eff_nom)
    delta_mu_MC = (muon_MC_eff_up - muon_MC_eff_nom)
    
    deff_MC = deff_MC_delectron_MC_eff_nom*delta_el_MC + deff_MC_dmuon_MC_eff_nom*delta_mu_MC
    
    deff_MC_rel = deff_MC/eff_MC 
    
    delta_SF = SF_nom*np.sqrt(deff_MC_rel**2 + deff_data_rel**2)

        # # up
    # eff_data_up  = Pass_mu_trig*muon_Data_eff_up + Pass_e_trig*electron_Data_eff_up - Pass_mu_trig*muon_Data_eff_up*Pass_e_trig*electron_Data_eff_up
    # eff_MC_up  = Pass_mu_trig*muon_MC_eff_up + Pass_e_trig*electron_MC_eff_up - Pass_mu_trig*muon_MC_eff_up*Pass_e_trig*electron_MC_eff_up
    # SF_up         =  eff_data_up/eff_MC_up
    # # down 
    # eff_data_down = Pass_mu_trig*muon_Data_eff_down + Pass_e_trig*electron_Data_eff_down - Pass_mu_trig*muon_Data_eff_down*Pass_e_trig*electron_Data_eff_down
    # eff_MC_down   = Pass_mu_trig*muon_MC_eff_down + Pass_e_trig*electron_MC_eff_down - Pass_mu_trig*muon_MC_eff_down*Pass_e_trig*electron_MC_eff_down
    # SF_down       =  eff_data_down/eff_MC_down

    #Saving the columns 
    events = set_ak_column_f32(events, f"Trigger_SF_nom", SF_nom)
    events = set_ak_column_f32(events, f"Trigger_SF_up", SF_nom + delta_SF )
    events = set_ak_column_f32(events, f"Trigger_SF_down", SF_nom - delta_SF)

    return events


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
