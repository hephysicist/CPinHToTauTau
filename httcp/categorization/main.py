# coding: utf-8

"""
Main categories file for the Higgs CP analysis
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import
from httcp.util import get_lep_p4, get_ip_p4

ak = maybe_import("awkward")
np = maybe_import("numpy")

#
# categorizer functions used by categories definitions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1

#Three general categories: etau, mutau and tautau
@categorizer(uses={'hcand_etau.*'})
def cat_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_etau.lep0.pt > 0, axis =1) > 0
    return events, mask 

@categorizer(uses={'hcand_mutau.*'})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_mutau.lep0.pt > 0, axis =1) > 0
    return events, mask

@categorizer(uses={'hcand_tautau.*'})
def cat_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_tautau.lep0.pt > 0, axis =1) > 0
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def os_charge(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
        mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].rel_charge < 0), axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def ss_charge(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
        mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].rel_charge > 0), axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def mt_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    mt_cut_value = self.config_inst.x.mt_cut_value
    for ch_str in channels:
        if ch_str != 'tautau':
            mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt <= mt_cut_value), axis=1),False)
            mask = mask & ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt >= 0), axis=1),False)
    #print(f"using mT cut {mt_cut_value}")
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def mt_inv_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    mt_cut_value = self.config_inst.x.mt_cut_value
    for ch_str in channels:
        if ch_str != 'tautau':
            mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt > mt_cut_value), axis=1),False)
            mask = mask & ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt < 200), axis=1),False)
    return events, mask

@categorizer(uses={"is_b_vetoed"})
def b_veto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ~events.is_b_vetoed #If event has b-veto is_b_vetoed mask is set to True, so to reject the event we need to inverse the mask
    return events, mask

@categorizer(uses={"is_b_vetoed"})
def b_veto_inv(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = events.is_b_vetoed #If event has b-veto is_b_vetoed mask is set to True, so to reject the event we need to inverse the mask
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def deep_tau_wp(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    
    deep_tau    =  self.config_inst.x.deep_tau
    wp_idx_ejet = deep_tau.vs_e_jet_wps
    wp_idx_mu   = deep_tau.vs_mu_wps
    
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for channel in channels:
        tau = events[f'hcand_{channel}'].lep1 
        channel_mask = ak.ones_like(events[f'hcand_{channel}'].lep1.rawIdx)
        if channel != 'tautau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= wp_idx_ejet[deep_tau.vs_jet[channel]])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= wp_idx_ejet[deep_tau.vs_e[channel]])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= wp_idx_mu[deep_tau.vs_mu[channel]])
        else:
            tau0 = events[f'hcand_{channel}'].lep0
            for the_tau in [tau, tau0]:
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSjet >= wp_idx_ejet[deep_tau.vs_jet[channel]])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSe   >= wp_idx_ejet[deep_tau.vs_jet[channel]])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSmu  >= wp_idx_mu[deep_tau.vs_mu[channel]])
        mask = mask | ak.fill_none(ak.firsts(channel_mask, axis=1),False)
    return events, mask


@categorizer(uses={'event', 'hcand_*'})
def deep_tau_wp_mtt(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    deep_tau_vs_e_jet_wps = self.config_inst.x.deep_tau.vs_e_jet_wps
    deep_tau_vs_mu_wps = self.config_inst.x.deep_tau.vs_mu_wps

    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for channel in channels:
        tau = events[f'hcand_{channel}'].lep1
        channel_mask = ak.ones_like(events[f'hcand_{channel}'].lep1.rawIdx)
        if channel == 'mutau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["Medium"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["Tight"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["Tight"])
        elif channel == 'etau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["Medium"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["Tight"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["Tight"])
        elif "tautau":
            tau0 = events[f'hcand_{channel}'].lep0
            for the_tau in [tau, tau0]:
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["Medium"])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["VVLoose"])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["VLoose"])
        mask = mask | ak.fill_none(ak.firsts(channel_mask, axis=1),False)
    return events, mask


@categorizer(uses={'event', 'hcand_*'})
def deep_tau_inv_wp(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    
    deep_tau    =  self.config_inst.x.deep_tau
    wp_idx_ejet = deep_tau.vs_e_jet_wps
    wp_idx_mu   = deep_tau.vs_mu_wps
    
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for channel in channels:
        tau = events[f'hcand_{channel}'].lep1 
        channel_mask = ak.ones_like(events[f'hcand_{channel}'].lep1.rawIdx)
        if channel != 'tautau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet < wp_idx_ejet[deep_tau.vs_jet[channel]])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= wp_idx_ejet['VVVLoose'])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= wp_idx_ejet[deep_tau.vs_e[channel]])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= wp_idx_mu[deep_tau.vs_mu[channel]])
        else:
            tau0 = events[f'hcand_{channel}'].lep0
            for the_tau in [tau, tau0]:
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSjet < wp_idx_ejet[deep_tau.vs_jet[channel]])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSjet >= wp_idx_ejet['VVVLoose'])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSe   >= wp_idx_ejet[deep_tau.vs_jet[channel]])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSmu  >= wp_idx_mu[deep_tau.vs_mu[channel]])
        mask = mask | ak.fill_none(ak.firsts(channel_mask, axis=1),False)
    return events, mask



@categorizer(uses={'event', 'hcand_*'})
def tau_endcap(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
            mask = mask | ak.fill_none(ak.firsts((np.abs(events[f'hcand_{ch_str}'].lep1.eta) > 1.2), axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def tau_barrel(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
            mask = mask | ak.fill_none(ak.firsts((np.abs(events[f'hcand_{ch_str}'].lep1.eta) <= 1.2), axis=1),False)
    return events, mask


@categorizer(uses={'event', 'hcand_*'})
def tau_eta2p3(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
            mask = mask | ak.fill_none(ak.firsts((np.abs(events[f'hcand_{ch_str}'].lep1.eta) < 2.3), axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def tau_no_fakes(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0] #We are processing a single channel at once
    if self.dataset_inst.is_mc:
        mask = ak.fill_none(ak.firsts(events[f'hcand_{channel}'].lep1.genPartFlav!=0, axis=1),False)
    else:
        mask = ak.ones_like(events.event, dtype=np.bool_)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def lep_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0] #We are processing a single channel at once
    if channel == 'etau': 
        isolation = events.hcand_etau.lep0.pfRelIso03_all < 0.3
    elif channel == 'mutau': 
        isolation = events.hcand_mutau.lep0.pfRelIso04_all < 0.15
    else:
        raise NotImplementedError(
                f'Can not find an isolation criteria for {channel} channel!')
    mask = ak.fill_none(ak.firsts(isolation, axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def lep_inv_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0] #We are processing a single channel at once
    if channel == 'etau': 
        isolation = events.hcand_etau.lep0.pfRelIso03_all >= 0.3
        upper_lim = events.hcand_etau.lep0.pfRelIso03_all < 2.
    elif channel == 'mutau': 
        isolation = events.hcand_mutau.lep0.pfRelIso04_all >= 0.15
        upper_lim = events.hcand_mutau.lep0.pfRelIso04_all < 0.5
    else:
        raise NotImplementedError(
                f'Can not find an isolation criteria for {channel} channel!')
    mask = ak.fill_none(ak.firsts(isolation, axis=1),False)
    mask = mask & ak.fill_none(ak.firsts(upper_lim, axis=1),False)
    return events, mask

from types import FunctionType
from copy import copy

def copy_function(fn, name):
    return FunctionType(
    copy(fn.__code__),
    copy(fn.__globals__),
    name=name,
    argdefs=copy(fn.__defaults__),
    closure=copy(fn.__closure__)
)

def single_pnet_dm(self: Categorizer, events: ak.Array, dm: int, **kwargs):
    channel = self.config_inst.channels.names()[0]
    mask = ak.fill_none(ak.firsts((events[f'hcand_{channel}'].lep1.decayModePNet == dm), axis=1),False)
    return events, mask

for the_dm in [0,1,2,10,11]:
    tmp_func = lambda self, events, **kwargs: single_pnet_dm(self, events, dm=the_dm, **kwargs)
    globals()[f'pnet_dm{the_dm}'] = categorizer(copy_function(tmp_func,f'pnet_dm{the_dm}'),
                                               uses={'event', 'hcand_*'}, )
    
def single_hps_dm(self: Categorizer, events: ak.Array, dm: int, **kwargs):
    channel = self.config_inst.channels.names()[0]
    mask = ak.fill_none(ak.firsts((events[f'hcand_{channel}'].lep1.decayMode == dm), axis=1),False)
    return events, mask

for the_dm in [0,1,10,11]:
    tmp_func = lambda self, events, **kwargs: single_hps_dm(self, events, dm=the_dm, **kwargs)
    globals()[f'hps_dm{the_dm}'] = categorizer(copy_function(tmp_func,f'hps_dm{the_dm}'),
                                               uses={'event', 'hcand_*'}, )

def _njets(self: Categorizer, events: ak.Array, njets, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = (events.n_jets==njets if njets < 2 else events.n_jets>=njets)
    return events, ak.fill_none(mask,False)

for the_njets in [0,1,2]:
    tmp_func = lambda self, events, **kwargs: _njets(self, events, njets=the_njets, **kwargs)
    globals()[f'njets_eq{the_njets}'] = categorizer(copy_function(tmp_func,f'njets_eq{the_njets}'),
                                               uses={'event', 'n_jets'}, )

#Helper category for check of jet veto maps (not used in the main analysis)
@categorizer(uses={"Jet*", "Muon*"}) 
def jet_veto_maps_jets(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    jet = events.Jet
    muon = events.Muon[events.Muon.isPFcand]
    jet_mask = (
        (jet.pt > 15) &
        (jet.jetId >= 2) &  # tight id 
        ((jet.chEmEF + jet.neEmEF) < 0.9) &
        ak.all(jet.metric_table(muon) >= 0.2, axis=2) &
        (np.abs(jet.eta) < 5.191)    # https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/JME_2022_Prompt_jetvetomaps.html
    )
    jet_mask = ak.fill_none(jet_mask, False)
    
    return events, ak.all(jet_mask, axis=1)

#Tau impact parameter cut used in definition of mupi category
@categorizer(uses={'event', 'hcand_*'})
def tau_ip_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0]
    mask = ak.fill_none(ak.firsts(events[f'hcand_{channel}'].lep1.ip_sig >= 1.25, axis=1),False)
    return events, mask

#Cut  on energy split between charged and neutral pion used in definition of murho and mu a1 categories
@categorizer(uses={'event', 'pion_E_split'})
def pion_E_split_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.fill_none(events.pion_E_split > 0.2, False)
    return events, mask

