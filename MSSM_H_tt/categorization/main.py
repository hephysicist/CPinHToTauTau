# coding: utf-8

"""
Main categories file for the Higgs CP analysis
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import


ak = maybe_import("awkward")
np = maybe_import("numpy")

#
# categorizer functions used by categories definitions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1

#Four general categories: etau, mutau, emu and tautau
@categorizer(uses={'hcand_etau.*'})
def cat_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_etau.lep0.pt > 0, axis =1) > 0
    return events, mask 

@categorizer(uses={'hcand_mutau.*'})
def cat_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_mutau.lep0.pt > 0, axis =1) > 0
    return events, mask

@categorizer(uses={'hcand_emu.*'})
def cat_emu(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = ak.num(events.hcand_emu.lep0.pt > 0, axis =1) > 0
    return events, mask 

@categorizer(uses={'event', 'hcand_*'})
def lep_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0] #We are processing a single channel at once
    if channel == 'emu': 
        isolation = events.hcand_emu.lep1.pfRelIso04_all < 0.2
    else:
        raise NotImplementedError(
                f'Can not find an isolation criteria for {channel} channel!')
    mask = ak.fill_none(ak.firsts(isolation, axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def lep_inv_iso(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channel = self.config_inst.channels.names()[0] #We are processing a single channel at once
    if channel == 'emu':
        isolation = events.hcand_emu.lep1.pfRelIso04_all >= 0.2
        upper_lim = events.hcand_emu.lep1.pfRelIso04_all < 0.5
    else:
        raise NotImplementedError(
                f'Can not find an isolation criteria for {channel} channel!')
    mask = ak.fill_none(ak.firsts(isolation, axis=1),False)
    mask = mask & ak.fill_none(ak.firsts(upper_lim, axis=1),False)
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
def mt_inv_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
        if ch_str != 'tautau':
            mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt > 50), axis=1),False)
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def mt_cut(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
        if ch_str != 'tautau':
            mask = mask | ak.fill_none(ak.firsts((events[f'hcand_{ch_str}'].mt <= 50), axis=1),False)
    return events, mask

@categorizer(uses={"N_b_jets"})
def Zero_b_jets(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = events.N_b_jets == 0 
    return events, mask

@categorizer(uses={"N_b_jets"})
def One_b_jets(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = events.N_b_jets == 1 
    return events, mask

@categorizer(uses={"N_b_jets"})
def At_least_2_b_jets(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = events.N_b_jets >= 2 
    return events, mask

@categorizer(uses={"OC_lepton_veto"})
def OC_lepton_veto(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = events.OC_lepton_veto
    return events, mask

@categorizer(uses={'event', 'hcand_*'})
def deep_tau_wp(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    channels = self.config_inst.channels.names()
    deep_tau_vs_e_jet_wps = self.config_inst.x.deep_tau.vs_e_jet_wps
    deep_tau_vs_mu_wps = self.config_inst.x.deep_tau.vs_mu_wps
    
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for channel in channels:
        tau = events[f'hcand_{channel}'].lep1 
        channel_mask = ak.ones_like(events[f'hcand_{channel}'].lep1.rawIdx)
        if channel == 'mutau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["Medium"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["VVLoose"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["Tight"])
        elif channel == 'etau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet >= deep_tau_vs_e_jet_wps["Medium"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["Tight"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["VLoose"])
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
    deep_tau_vs_e_jet_wps = self.config_inst.x.deep_tau.vs_e_jet_wps
    deep_tau_vs_mu_wps = self.config_inst.x.deep_tau.vs_mu_wps
    
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for channel in channels:
        tau = events[f'hcand_{channel}'].lep1 
        channel_mask = ak.ones_like(events[f'hcand_{channel}'].lep1.rawIdx)
        if channel == 'mutau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet < deep_tau_vs_e_jet_wps["Medium"]) #This cut is reversed
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["VVLoose"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["Tight"])
        elif channel == 'etau':
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSjet < deep_tau_vs_e_jet_wps["Medium"]) #This cut is reversed
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["Tight"])
            channel_mask = channel_mask & (tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["VLoose"])
        elif "tautau":
            tau0 = events[f'hcand_{channel}'].lep0
            for the_tau in [tau, tau0]:
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSjet < deep_tau_vs_e_jet_wps["Medium"]) #This cut is reversed
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSe   >= deep_tau_vs_e_jet_wps["VVLoose"])
                channel_mask = channel_mask & (the_tau.idDeepTau2018v2p5VSmu  >= deep_tau_vs_mu_wps["VLoose"])
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

@categorizer(uses={'D_zeta'})
def D_zeta_cut_low(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = (events.D_zeta >= -35) & (events.D_zeta < -10)
    return events, mask

@categorizer(uses={'D_zeta'})
def D_zeta_cut_mid(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = (events.D_zeta >= -10) & (events.D_zeta < 30)
    return events, mask

@categorizer(uses={'D_zeta'})
def D_zeta_cut_high(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    mask = (events.D_zeta >= 30)
    return events, mask
