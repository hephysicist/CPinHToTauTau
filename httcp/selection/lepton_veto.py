# coding: utf-8

"""
Extra-Lepton-Veto
https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from law.util import DotDict



np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass",
    },
    exposed=False,
)
def extra_lepton_veto(
        self: Selector,
        events: ak.Array,
        extra_electron_index: ak.Array,
        extra_muon_index: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    #Selectin ALL leptons that pass extra lepton kinematic cuts
    
    extra_lep  = ak.Array(ak.concatenate([events.Muon[extra_muon_index],
                                          events.Electron[extra_electron_index]], axis=-1),
                          behavior=coffea.nanoevents.methods.nanoaod.behavior)
    extra_lep  = ak.with_name(extra_lep, "PtEtaPhiMLorentzVector")
    
    #check if the events contain a single pair i.e. number of rawIdxs in the hcand object is equal to 2
    # has_single_pair = ak.sum(ak.num(hcand_pair.rawIdx, axis=2),axis=1) == 2 ATTENTION: hcand has a single pair by construction 
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    has_no_exra_lep = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels: 
        hcand = events[f'hcand_{ch_str}']
        for lep_str in hcand.fields:
            lep = hcand[lep_str]
            #Calculate dR between first lepton and electrons or muons from loose selection list
            #Attention: this list also contains dR between lep1 and itself, so all entries have at least one dR that is 0 or a tiny number
            delta_r = ak.flatten(extra_lep.metric_table(lep),axis=2)
            n_sep_leps = ak.sum(delta_r > 0.5, axis=1)
            if ch_objects[ch_str][lep_str] == 'Tau':
                 #For tau objects there also exists a lepton from the pair, so it will trigger delta_r cut, so to tirgger the veto there should be at least two particles
                lep_mask = n_sep_leps > 1
            else:
                #For electrons and muons there we need to check if there exists a lepton, separated from the pair lepton
                lep_mask = n_sep_leps > 0
            has_no_exra_lep = has_no_exra_lep | lep_mask
    return events, SelectionResult(
        steps={"extra_lepton_veto": ~has_no_exra_lep #Since it's veto, this field shoud be reversed 
               })



@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass", "Muon.charge",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass", "Electron.charge",
    },
    exposed=False,
)
def double_lepton_veto(
        self: Selector,
        events: ak.Array,
        double_veto_electron_index: ak.Array,
        double_veto_muon_index: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    double_veto_muon     = ak.Array(events.Muon[double_veto_muon_index], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)

    double_veto_muon     = ak.with_name(double_veto_muon, "PtEtaPhiMLorentzVector")

    double_veto_electron = ak.Array(events.Electron[double_veto_electron_index], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)
    double_veto_electron = ak.with_name(double_veto_electron, "PtEtaPhiMLorentzVector")

    mu_pair = ak.combinations(double_veto_muon, 2, axis=1)
    el_pair = ak.combinations(double_veto_electron, 2, axis=1)

    mu1,mu2 = ak.unzip(mu_pair)
    el1,el2 = ak.unzip(el_pair)

    presel_mask = lambda leps1, leps2: ((leps1.charge * leps2.charge < 0) & (leps1.delta_r(leps2) > 0.15))

    dlveto_mu_mask = presel_mask(mu1,mu2)
    dlveto_el_mask = presel_mask(el1,el2)

    dl_mu_veto = ak.sum(dlveto_mu_mask, axis=1) == 0
    dl_el_veto = ak.sum(dlveto_el_mask, axis=1) == 0

    dl_veto    = dl_mu_veto & dl_el_veto

    return events, SelectionResult(steps={"dilepton_veto": dl_veto})
    
@selector(
    uses={
        "LHEPart.pdgId",
    },
    exposed=False,
)
def tau_veto(
        self: Selector,
        events: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    n_taus = ak.sum(abs(events.LHEPart.pdgId) == 15, axis=-1)
    return events, SelectionResult(steps={"tau_veto": (n_taus == 0)})