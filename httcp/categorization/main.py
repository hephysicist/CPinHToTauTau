# coding: utf-8

"""
Main categories file for the Higgs CP analysis
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


#
# categorizer functions used by categories definitions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1


@categorizer(uses={"Jet.pt"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2


# ---------------------------------------------------------- #
#                            e-tau                           #
# ---------------------------------------------------------- #
@categorizer(uses={"channel_id"})
def sel_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    return events, events["channel_id"] == ch.id

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_etau_pion(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 0), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_etau_rho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 1), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_etau_a1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("etau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 10), axis=1) == 1
    return events, ch_mask & dm_mask

# ---------------------------------------------------------- #
#                             mu-tau                         #
# ---------------------------------------------------------- #

@categorizer(uses={"channel_id"})
def sel_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    return events, events["channel_id"] == ch.id

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_mutau_pion(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 0), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_mutau_rho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 1), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_mutau_a1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("mutau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode[:,1:2] == 10), axis=1) == 1
    return events, ch_mask & dm_mask


# ---------------------------------------------------------- #
#                            tau-tau                         #
# ---------------------------------------------------------- #

@categorizer(uses={"channel_id"})
def sel_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    return events, events["channel_id"] == ch.id

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_pionpion(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode == 0), axis=1) == 2
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_rhorho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum((events.hcand.decayMode == 1), axis=1) == 2
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_a1a1(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    dm = events.hcand.decayMode == 10
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum(dm, axis=1) == 2
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_pionrho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    dm_hcand1 = events.hcand.decayMode[:,0:1]
    dm_hcand2 = events.hcand.decayMode[:,1:2]    
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum(((dm_hcand1 == 0) & (dm_hcand2 == 1)) | ((dm_hcand1 == 1) & (dm_hcand2 == 0)), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_a1pion(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    dm_hcand1 = events.hcand.decayMode[:,0:1]
    dm_hcand2 = events.hcand.decayMode[:,1:2]    
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum(((dm_hcand1 == 0) & (dm_hcand2 == 10)) | ((dm_hcand1 == 10) & (dm_hcand2 == 0)), axis=1) == 1
    return events, ch_mask & dm_mask

@categorizer(uses={"channel_id", "hcand.decayMode"})
def sel_tautau_a1rho(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    ch = self.config_inst.get_channel("tautau")
    dm_hcand1 = events.hcand.decayMode[:,0:1]
    dm_hcand2 = events.hcand.decayMode[:,1:2]
    ch_mask = events["channel_id"] == ch.id
    dm_mask = ak.sum(((dm_hcand1 == 10) & (dm_hcand2 == 1)) | ((dm_hcand1 == 1) & (dm_hcand2 == 10)), axis=1) == 1
    return events, ch_mask & dm_mask



##################################
### Cathegories for ABCD method### 
##################################


@categorizer(uses={"hcand.rel_charge", "Muon.pfRelIso04_all"})
def cat_c(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Control region ( iso < 0.15, same sign pair)
    sel = (events.hcand.rel_charge[:,0] > 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) < 0.15)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_d(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Signal region ( iso < 0.15, opposite sign pair)
    sel = (events.hcand.rel_charge[:,0] < 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) < 0.15)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_a(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Region for transfer factor calculation( iso > 0.15, same sign pair)
    sel = (events.hcand.rel_charge[:,0] > 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) >= 0.15)  \
        & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) <= 0.30)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Region for transfer factor calculation( iso > 0.15, opposite sign pair)
    sel = (events.hcand.rel_charge[:,0] < 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) >= 0.15) \
        & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) <= 0.30)
    return events, sel
