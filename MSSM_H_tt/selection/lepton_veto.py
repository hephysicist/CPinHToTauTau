# coding: utf-8

"""
Extra-Lepton-Veto
https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2325&ancode=HIG-20-006&tp=an&line=HIG-20-006
http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_192_v15.pdf
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
 
@selector(
    exposed=False,
)
def single_lepton_veto(
        self: Selector,
        events: ak.Array,
        single_electron_index: ak.Array,
        single_muon_index: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    channels = self.config_inst.channels.names()
    number_of_muons =  ak.num(single_muon_index)
    number_of_electrons =  ak.num(single_electron_index)
    
    # Use this in etau channel 
    single_muon_veto = number_of_muons == 0
    # Use this in mutau channel
    single_electron_veto = number_of_electrons == 0   
    # Use this in tautau channel
    single_ele_mu_veto = single_electron_veto & single_muon_veto
    for ch_str in channels:
        if ch_str == 'etau':
            single_lepton_veto = single_muon_veto
        elif ch_str == 'mutau':
            single_lepton_veto = single_electron_veto
        if ch_str == 'emu':
            single_lepton_veto = number_of_electrons >=0 #We don't have this veto in emu, so it's always True
        elif ch_str == 'tautau':
            single_lepton_veto = single_ele_mu_veto
    return events, SelectionResult(
        steps={"single_lepton_veto": single_lepton_veto})
    
@selector(
    exposed=False,
)
def second_lepton_veto(
        self: Selector,
        events: ak.Array,
        single_electron_index: ak.Array,
        single_muon_index: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    channels = self.config_inst.channels.names()
    number_of_muons =  ak.num(single_muon_index)
    number_of_electrons =  ak.num(single_electron_index)
    
    # Use this in mutau channel 
    second_muon_veto = number_of_muons >= 2
    # Use this in etau channel
    second_electron_veto = number_of_electrons >= 2  
    # Use this in tautau channel, always True because we don't apply this veto
    second_ele_mu_veto = number_of_muons >= 0
    for ch_str in channels:
        if ch_str == 'etau':
            second_lepton_veto = ~second_electron_veto
        elif ch_str == 'mutau':
            second_lepton_veto = ~second_muon_veto
        if ch_str == 'emu':
            second_lepton_veto = ((~second_electron_veto) & (~second_muon_veto))
        elif ch_str == 'tautau':
            second_lepton_veto = second_ele_mu_veto
    return events, SelectionResult(
        steps={"second_lepton_veto": second_lepton_veto})

@selector(
    uses={
        "Muon.pt", "Muon.eta", "Muon.phi", "Muon.mass", "Muon.charge",
        "Electron.pt", "Electron.eta", "Electron.phi", "Electron.mass", "Electron.charge",
    },
    exposed=False,
)
def OC_lepton_veto(
        self: Selector,
        events: ak.Array,
        OC_veto_electron_indices: ak.Array,
        OC_veto_muon_indices: ak.Array,
        **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    
    channels = self.config_inst.channels.names()
        
    OC_veto_electron = ak.Array(events.Electron[OC_veto_electron_indices], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)
    OC_veto_electron = ak.with_name(OC_veto_electron, "PtEtaPhiMLorentzVector")
    
    OC_veto_muon     = ak.Array(events.Muon[OC_veto_muon_indices], 
                                    behavior=coffea.nanoevents.methods.nanoaod.behavior)
    OC_veto_muon     = ak.with_name(OC_veto_muon, "PtEtaPhiMLorentzVector")
   
    el_pair = ak.combinations(OC_veto_electron, 2, axis=1)
    mu_pair = ak.combinations(OC_veto_muon, 2, axis=1)
   
    el1,el2 = ak.unzip(el_pair)
    mu1,mu2 = ak.unzip(mu_pair)

    presel_mask = lambda leps1, leps2: ((leps1.charge * leps2.charge < 0) & (leps1.delta_r(leps2) > 0.15))
    
    OC_veto_el_mask = presel_mask(el1,el2)
    OC_veto_mu_mask = presel_mask(mu1,mu2)
    
    # Use this in etau channel
    OC_veto_ele  = ak.sum(OC_veto_el_mask, axis=1) == 0
    # Use this in mutau channel 
    OC_veto_mu   =  ak.sum(OC_veto_mu_mask, axis=1) == 0
    # Use this in emu and tautau channel, always True because we don't apply this veto
    OC_ele_mu_veto = ak.num(OC_veto_electron_indices) >= 0
    
    for ch_str in channels:
        if ch_str == 'etau':
            OC_lepton_veto = OC_veto_ele
        elif ch_str == 'mutau':
            OC_lepton_veto = OC_veto_mu
        if ch_str == 'emu':
            OC_lepton_veto = OC_ele_mu_veto
        elif ch_str == 'tautau':
            OC_lepton_veto = OC_ele_mu_veto
    events = set_ak_column(events, "OC_lepton_veto", OC_lepton_veto) 
    return events, SelectionResult(
        steps={"OC_lepton_veto": OC_lepton_veto})