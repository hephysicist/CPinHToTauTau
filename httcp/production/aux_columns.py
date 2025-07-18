"""
Produce channel_id column. This function is called in the main selector
"""
import functools
import law
from columnflow.production import Producer, producer
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT
from columnflow.util import maybe_import
from law.util import DotDict
from httcp.util import get_lep_p4, get_vec_p3, to_pt_eta_phi_m

np = maybe_import("numpy")
ak = maybe_import("awkward")

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

logger = law.logger.get_logger(__name__)


@producer(
    produces={
        "channel_id",
    },
    exposed=False,
)
def channel_id(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    # Create a template array filled with zeros
    channel_id = ak.zeros_like(ak.local_index(events.event), dtype=np.uint8)

    # Create a mask array to check for channel orthogonality
    and_mask = ak.ones_like(ak.local_index(events.event), dtype=np.bool_)

    for channel in self.config_inst.channels.names():
        the_mask = ak.num(events[f'hcand_{channel}'].lep0, axis=1) > 0
        and_mask = and_mask & the_mask
        channel_id = ak.where(the_mask,
                              self.config_inst.get_channel(channel).id,
                              channel_id)

    channel_id = ak.values_astype(channel_id, np.uint8)
    events = set_ak_column(events, "channel_id", channel_id)

    return events


@producer(
    uses={"Jet.*"},
    produces={"Jet.*"},
    exposed=False,
)
def create_jetID_masks(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    For nanoaod v13 and v14 jetID is bugged, so there is a special procedure to apply Tight jet ID 
    """
    nano_version = self.config_inst.campaign.x.version
    jets = events.Jet
    if nano_version in [13, 14]:
        print(f'Applying custom tightJetID for nanoAOD v{nano_version}...')
        tightID_eta_2p6 = ((jets.neHEF < 0.99)
                           & (jets.neEmEF < 0.9)
                           & ((jets.chMultiplicity + jets.neMultiplicity) > 1)
                           & (jets.chHEF > 0.01)
                           # Tight criteria for |eta| < 2.6
                           & (jets.chMultiplicity > 0))

        # Tight criteria for 2.6 < |eta| <= 2.7
        tightID_eta_2p6_to_2p7 = ((jets.neHEF < 0.9) & (jets.neEmEF < 0.99))

        # Tight criteria for 2.7 < |eta| <= 3.0
        tightID_eta_2p7_to_3p0 = (jets.neHEF < 0.99)

        tightID_eta_geq_3p0 = ((jets.neMultiplicity >= 2) & (
            jets.neEmEF < 0.4))  # Tight criteria for |eta| >= 3.0

        pass_tightID = (((abs(jets.eta) < 2.6) & tightID_eta_2p6)
                        | ((abs(jets.eta) >= 2.6) & (abs(jets.eta) < 2.7) & tightID_eta_2p6_to_2p7)
                        | ((abs(jets.eta) >= 2.7) & (abs(jets.eta) < 3.0) & tightID_eta_2p7_to_3p0)
                        | ((abs(jets.eta) >= 3.0) & tightID_eta_geq_3p0)
                        )

        pass_tightID_lep_veto = ak.where((abs(jets.eta) < 2.7),
                                         pass_tightID & (jets.muEF < 0.8) & (
                                             jets.chEmEF < 0.8),
                                         pass_tightID)
    else:
        pass_tightID_lep_veto = (abs(jets.eta) < 2.7) & ((jets.jetId & 3) != 0)

    jets["pass_tightID"] = pass_tightID
    jets["pass_tightID_lep_veto"] = pass_tightID_lep_veto
    events = set_ak_column(events, "Jet", jets)
    return events


@producer(
    uses={f"Jet.{var}" for var in
          [
        "pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"
    ]} | {"hcand_*"},
    produces={"is_b_vetoed"},
    exposed=False,
)
def jet_veto(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces 'is_b_vetoed' column to be used in categorisation
    """
    year = self.config_inst.x.year
    tag = self.config_inst.x.tag
    btag_wp = self.config_inst.x.btag_working_points[year][tag].deepjet.medium

    # nominal jet selection
    jet_pt_sorted_idx = ak.argsort(events.Jet.pt, axis=1, ascending=False)
    sorted_jets = events.Jet[jet_pt_sorted_idx]
    jet_selections = {
        "jet_pt_20": sorted_jets.pt > 20.0,
        "jet_eta_2.5": abs(sorted_jets.eta) < 2.5,
        "jet_id": sorted_jets.pass_tightID_lep_veto,
        "btag_wp_medium": sorted_jets.btagDeepFlavB >= btag_wp
    }
    jet_obj_mask = ak.ones_like(jet_pt_sorted_idx, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel
    event_mask = ak.zeros_like(events.event, dtype=np.bool_)
    for the_ch in self.config_inst.channels.names():
        hcand = events[f'hcand_{the_ch}']
        for lep_str in [field for field in hcand.fields if 'lep' in field]:
            lep = ak.firsts(hcand[lep_str])
            seed_idx = ak.fill_none(lep.jetIdx, -1)
            jet_obj_mask = jet_obj_mask & (jet_pt_sorted_idx != seed_idx)
            presel_jet = ak.drop_none(ak.mask(sorted_jets, jet_obj_mask))
            jet_tau_pairs = ak.cartesian([presel_jet, lep], axis=1)
            jet_br, lep_br = ak.unzip(jet_tau_pairs)
            delta_r = jet_br.delta_r(lep_br)
            # make joint mask for first and second tau.
            event_mask = event_mask | ak.any(delta_r < 0.4, axis=1)
    events = set_ak_column(events, "is_b_vetoed", event_mask)
    return events


@producer(
    uses={f"Jet.{var}" for var in
          ["pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"
           ]} | {"hcand_*"},
    produces={
        "n_jets",
        "lead_jet.*",
        "sublead_jet.*",
        "dijet.*",
    },
    exposed=False,
)
def jet_pt_def(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces a mask to plot the jet pt observable
    """
    # nominal jet selection
    jet_pt_sorted_idx = ak.argsort(events.Jet.pt, axis=1, ascending=False)
    sorted_jets = events.Jet[jet_pt_sorted_idx]
    jet_selections = {
        "jet_pt_20": sorted_jets.pt > 20.0,
        "jet_eta_4.7": abs(sorted_jets.eta) < 4.7,
        "jet_id": sorted_jets.pass_tightID_lep_veto,
    }
    jet_obj_mask = ak.ones_like(jet_pt_sorted_idx, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel

    ch_str = self.config_inst.channels.names()[0]
    hcand = events[f'hcand_{ch_str}']
    mask = ak.ones_like(sorted_jets.pt, dtype=np.bool_)
    empty_p4 = ak.zeros_like(get_lep_p4(hcand.lep0))
    for lep_str in ['lep0', 'lep1']:
        lep = ak.firsts(hcand[lep_str])
        seed_idx = ak.fill_none(lep.jetIdx, -1)
        jet_obj_mask_seed_idx = jet_obj_mask & (jet_pt_sorted_idx != seed_idx)
        presel_jet = ak.mask(sorted_jets, jet_obj_mask_seed_idx)
        _jet_br, _lep_br = ak.unzip(ak.cartesian([presel_jet, lep], axis=1))
        jet_p4 = get_lep_p4(_jet_br)
        lep_p4 = get_lep_p4(_lep_br)
        mask = mask & ak.fill_none((jet_p4.delta_r(lep_p4) > 0.4), False)

    jets = get_lep_p4(sorted_jets[mask])
    jet_pt_mask = ((jets.pt > 20) & (abs(jets.eta) <= 2.5))
    jet_pt_mask = jet_pt_mask | ((jets.pt > 50) & (
        (abs(jets.eta) <= 3.0) & (abs(jets.eta) > 2.5)))
    jet_pt_mask = jet_pt_mask | ((jets.pt > 30) & ((abs(jets.eta) > 3.0)))

    sel_jets = ak.drop_none(ak.mask(jets, jet_pt_mask))
    njets = ak.num(sel_jets, axis=1)

    lead_jet_p4 = ak.where(njets > 0,
                           sel_jets[:, :1],
                           empty_p4)

    sublead_jet_p4 = ak.where(njets > 1,
                              sel_jets[:, 1:2],
                              empty_p4)

    dijet_p4 = to_pt_eta_phi_m(lead_jet_p4 + sublead_jet_p4)

    lead_jet = {}
    lead_jet_mask = (lead_jet_p4.pt > 20)
    sublead_jet = {}
    sublead_jet_mask = (sublead_jet_p4.pt > 20)
    dijet = {}
    dijet_mask = (dijet_p4.pt > 20)
    empty_float = lambda arr : ak.full_like(arr, EMPTY_FLOAT)

    for var in ['pt', 'eta', 'phi', 'mass']: 
        lead_jet[var] = ak.where(lead_jet_mask,
                                 getattr(lead_jet_p4, var),
                                 empty_float(getattr(lead_jet_p4, var))
                                 )
        sublead_jet[var] = ak.where(sublead_jet_mask,
                                    getattr(sublead_jet_p4, var),
                                    empty_float(getattr(sublead_jet_p4, var))
                                    )
        dijet[var] = ak.where(dijet_mask,
                              getattr(dijet_p4, var),
                              empty_float(getattr(dijet_p4, var))
                              )

    for func_name in ['deltaeta', 'deltaphi', 'delta_r']:
        func = getattr(lead_jet_p4, func_name)
        dijet[func_name] = ak.where(njets >= 2,
                                    func(sublead_jet_p4),
                                    empty_float(dijet_p4.pt)
                                    )

    events = set_ak_column(events, 'lead_jet', ak.zip(lead_jet))
    events = set_ak_column(events, 'sublead_jet', ak.zip(sublead_jet))
    events = set_ak_column(events, 'dijet', ak.zip(dijet))
    events = set_ak_column(events, 'n_jets', njets)
    return events


@producer(
    uses={f"Jet.{var}" for var in
          ["pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"
           ]} | {"hcand_*"},
    produces={"n_jets_tag"},
    exposed=False,
)
def jets_taggable(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces a mask to plot the jet pt observable
    """
    # nominal jet selection
    jet_pt_sorted_idx = ak.argsort(events.Jet.pt, axis=1, ascending=False)
    sorted_jets = events.Jet[jet_pt_sorted_idx]
    jet_selections = {
        "jet_pt_20": sorted_jets.pt > 20.0,
        "jet_eta_2.5": abs(sorted_jets.eta) < 2.5,
        "jet_id": sorted_jets.pass_tightID_lep_veto,
    }
    jet_obj_mask = ak.ones_like(jet_pt_sorted_idx, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel

    for the_ch in self.config_inst.channels.names():
        hcand = events[f'hcand_{the_ch}']

        for lep_str in [field for field in hcand.fields if 'lep' in field]:
            lep = ak.firsts(hcand[lep_str])
            seed_idx = ak.fill_none(lep.jetIdx, -1)
            jet_obj_mask_seed_idx = jet_obj_mask & (
                jet_pt_sorted_idx != seed_idx)
            presel_jet = ak.drop_none(
                ak.mask(sorted_jets, jet_obj_mask_seed_idx))
            jet_tau_pairs = ak.cartesian([presel_jet, lep], axis=1)
            jet_br, lep_br = ak.unzip(jet_tau_pairs)
            delta_phi = (jet_br.phi - lep_br.phi)
            delta_phi = ak.where(
                delta_phi > np.pi, delta_phi - 2*np.pi, delta_phi)
            delta_phi = ak.where(delta_phi < -np.pi,
                                 delta_phi + 2*np.pi, delta_phi)
            delta_eta = (jet_br.eta - lep_br.eta)
            delta_r = np.sqrt(delta_phi**2 + delta_eta**2)
            jet_br_to_plot = jet_br[delta_r > 0.4]
            if lep_str == 'lep0':
                jet_vs_lep0 = jet_br[delta_r > 0.4]
            elif lep_str == 'lep1':
                jet_vs_lep1 = jet_br[delta_r > 0.4]

    lep0_max_obj = ak.max(ak.num(jet_vs_lep0.pt))
    lep1_max_obj = ak.max(ak.num(jet_vs_lep1.pt))
    max_len = ak.max([lep0_max_obj, lep1_max_obj])

    jet_vs_lep0_pt = ak.pad_none(jet_vs_lep0.pt, max_len)
    jet_vs_lep0_pt_tag = ak.fill_none(jet_vs_lep0_pt, -999)
    jet_vs_lep1_pt = ak.pad_none(jet_vs_lep1.pt, max_len)
    jet_vs_lep1_pt_tag = ak.fill_none(jet_vs_lep1_pt, -999)

    n_jets_vs_lep0_tag = ak.sum(jet_vs_lep0_pt_tag > 20, axis=1)
    n_jets_vs_lep1_tag = ak.sum(jet_vs_lep1_pt_tag > 20, axis=1)
    n_jets_mask_tag = (n_jets_vs_lep0_tag == n_jets_vs_lep1_tag)
    n_jets_taggable = ak.where(n_jets_mask_tag, n_jets_vs_lep0_tag, 0)
    events = set_ak_column(events, "n_jets_tag", n_jets_taggable)
    return events


@producer(
    uses={f"Jet.{var}" for var in
          [
        "pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"
    ]} | {"hcand_*"},
    produces={"N_b_jets"},
    exposed=False,
)
def number_b_jet(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces 'N_b_jets' column to be used in categorisation
    """
    year = self.config_inst.x.year
    tag = self.config_inst.x.tag
    btag_wp = self.config_inst.x.btag_working_points[year][tag].deepjet.medium
    # nominal jet selection
    jet_pt_sorted_idx = ak.argsort(events.Jet.pt, axis=1, ascending=False)
    sorted_jets = events.Jet[jet_pt_sorted_idx]
    jet_selections = {
        "jet_pt_20": sorted_jets.pt > 20.0,
        "jet_eta_2.5": abs(sorted_jets.eta) < 2.5,
        "jet_id": sorted_jets.pass_tightID_lep_veto,
        "btag_wp_medium": sorted_jets.btagDeepFlavB >= btag_wp
    }
    jet_obj_mask = ak.ones_like(jet_pt_sorted_idx, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel

    for the_ch in self.config_inst.channels.names():
        hcand = events[f'hcand_{the_ch}']

        for lep_str in [field for field in hcand.fields if 'lep' in field]:
            lep = ak.firsts(hcand[lep_str])
            seed_idx = ak.fill_none(lep.jetIdx, -1)
            jet_obj_mask_seed_idx = jet_obj_mask & (
                jet_pt_sorted_idx != seed_idx)
            presel_jet = ak.drop_none(
                ak.mask(sorted_jets, jet_obj_mask_seed_idx))
            jet_tau_pairs = ak.cartesian([presel_jet, lep], axis=1)
            jet_br, lep_br = ak.unzip(jet_tau_pairs)
            delta_phi = jet_br.phi - lep_br.phi
            delta_phi = ak.where(
                delta_phi > np.pi, delta_phi - 2*np.pi, delta_phi)
            delta_phi = ak.where(delta_phi < -np.pi,
                                 delta_phi + 2*np.pi, delta_phi)
            delta_eta = jet_br.eta - lep_br.eta
            delta_r = np.sqrt(delta_phi**2 + delta_eta**2)
            jet_br_to_plot = jet_br[delta_r > 0.4]
            if lep_str == 'lep0':
                jet_vs_lep0 = jet_br[delta_r > 0.4]
            elif lep_str == 'lep1':
                jet_vs_lep1 = jet_br[delta_r > 0.4]

    lep0_max_obj = ak.max(ak.num(jet_vs_lep0.pt))
    lep1_max_obj = ak.max(ak.num(jet_vs_lep1.pt))
    max_len = ak.max([lep0_max_obj, lep1_max_obj])

    jet_vs_lep0_pt = ak.pad_none(jet_vs_lep0.pt, max_len)
    jet_vs_lep0_pt_tag = ak.fill_none(jet_vs_lep0_pt, -999)
    jet_vs_lep1_pt = ak.pad_none(jet_vs_lep1.pt, max_len)
    jet_vs_lep1_pt_tag = ak.fill_none(jet_vs_lep1_pt, -999)

    n_jets_vs_lep0_tag = ak.sum((jet_vs_lep0_pt_tag > 20), axis=1)
    n_jets_vs_lep1_tag = ak.sum((jet_vs_lep1_pt_tag > 20), axis=1)
    n_jets_mask_tag = (n_jets_vs_lep0_tag == n_jets_vs_lep1_tag)
    n_jets_taggable = ak.where(n_jets_mask_tag, n_jets_vs_lep0_tag, 0)

    events = set_ak_column(events, "N_b_jets", n_jets_taggable)
    return events


def is_pion(prods): return ((np.abs(prods.pdgId) == 211))


def is_photon(prods): return prods.pdgId == 22


def has_one_pion(prods): return (ak.sum(is_pion(prods),   axis=1) == 1)


def has_three_pions(prods): return (ak.sum(is_pion(prods),   axis=1) == 3)


def has_photons(prods): return (ak.sum(is_photon(prods), axis=1) > 0)


def has_no_photons(prods): return (ak.sum(is_photon(prods), axis=1) == 0)


@producer(
    uses={
        "TauProd.pdgId",
    },
    produces={
        "TauProd.mass", "TauProd.charge",
    },
    exposed=False,
)
def assign_tauprod_mass_charge(
        self: Selector,
        events: ak.Array,
        **kwargs
) -> ak.Array:

    # https://pdg.lbl.gov/2023/listings/particle_properties.html

    part_dict = DotDict.wrap({
        "pion_pm": {'mass': 0.13957039,  # GeV
                    'pdg_id': 211},
        "pion_0": {'mass': 0.1349768,  # GeV
                   'pdg_id': 111},
        "gamma": {'mass': 0.0,  # GeV
                  'pdg_id': 22},
        "kaon_pm": {'mass': 0.493677,  # GeV
                    'pdg_id': 321},
        "ele_pm": {'mass': 0.00051099895,  # GeV
                   'pdg_id': 11},
        "muon_pm": {'mass': 0.1056583755,  # GeV
                    'pdg_id': 13},
    })
    mass = ak.zeros_like(ak.local_index(
        events.TauProd.pdgId), dtype=np.float32)
    charge = ak.zeros_like(ak.local_index(
        events.TauProd.pdgId), dtype=np.int32)
    for part in part_dict:

        prod_id = events.TauProd.pdgId
        mass = ak.where(np.abs(prod_id) == part_dict[part].pdg_id,
                        part_dict[part].mass,
                        mass)
        if '_pm' in part:
            charge = ak.where(np.abs(prod_id) == part_dict[part].pdg_id,
                              np.sign(prod_id),
                              charge)

    events = set_ak_column_f32(events, "TauProd.mass", mass)
    events = set_ak_column_i32(events, "TauProd.charge", charge)
    return events


@producer(
    uses={
        "TauProd.*", assign_tauprod_mass_charge,
    },
    produces={
        "tau_decay_prods_*",
        assign_tauprod_mass_charge,
    },
    exposed=False,
)
def add_tau_prods(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, SelectionResult]:

    events = self[assign_tauprod_mass_charge](events)
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    mask = ak.zeros_like(events.event, dtype=np.bool_)
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        tau_decay_prods_dict = {}
        for lep_str in hcand.fields:
            if ch_objects[ch_str][lep_str] == 'Tau':
                tau = hcand[lep_str]
                tauprods = events.TauProd
                tau_idx = tau.rawIdx
                tauprods_mask = tauprods.tauIdx
                idx_pairs = ak.cartesian(
                    [tau_idx, tauprods.tauIdx], axis=1, nested=True)
                tau_idx_b, tau_prod_idx_b = ak.unzip(idx_pairs)
                tau2prod_match_mask = ak.fill_none(
                    tau_idx_b == tau_prod_idx_b, False)
                matched_tau_prods = tauprods[ak.flatten(
                    tau2prod_match_mask, axis=2)]
                # DM0
                tau = ak.firsts(tau)
                mask = mask | ak.fill_none(
                    tau.decayMode == 0, False) & has_one_pion(matched_tau_prods)#TODO: release this mask, make n_pions >=1
                # DM1
                mask = mask | ak.fill_none(tau.decayMode == 1, False) & has_one_pion(
                    matched_tau_prods) & has_photons(matched_tau_prods)#TODO: release this mask, make n_pions >=1
                # DM10
                mask = mask | ak.fill_none(
                    tau.decayMode == 10, False) & has_three_pions(matched_tau_prods)#TODO: release this mask, make n_pions >=3
                # DM11
                mask = mask | ak.fill_none(tau.decayMode == 11, False) & has_three_pions(
                    matched_tau_prods) & has_photons(matched_tau_prods)#TODO: release this mask, make n_pions >=3
                events = set_ak_column(
                    events, f'tau_decay_prods_{ch_str}_{lep_str}',  matched_tau_prods)
            else:
                pass

    return events, SelectionResult(
        steps={
            "decay_prods_are_ok": mask,
        },
    )


def egamma_mask(tauprod): return (
    (np.abs(tauprod.pdgId) == 11) + (tauprod.pdgId == 22))


def pion_mask(tauprod): return np.abs(tauprod.pdgId) == 211


@producer(
    uses={
        "tau_decay_prods_*",
    },
    produces={
        "pion_E_split",
    },
    exposed=False,
)
def pion_energy_split(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    channel = self.config_inst.channels.names()[0]
    tauprods = events[f'tau_decay_prods_{channel}_lep1']
    charged_pion_mask = pion_mask(tauprods)
    em_mask = egamma_mask(tauprods)
    charged_pion = ak.firsts(get_lep_p4(tauprods[charged_pion_mask]), axis=1)
    neutral_pion = get_lep_p4(tauprods[em_mask]).sum()

    mask = (ak.num(charged_pion_mask, axis=1) > 0) & (
        ak.num(em_mask, axis=1) > 0)
    mask = mask & (charged_pion.E > 0) & (neutral_pion.E > 0)

    pion_E_split = ak.where(mask,
                            np.abs(charged_pion.E - neutral_pion.E) /
                            (charged_pion.E + neutral_pion.E),
                            EMPTY_FLOAT)
    pion_E_split = ak.fill_none(pion_E_split, EMPTY_FLOAT)
    events = set_ak_column_f32(events, "pion_E_split", pion_E_split)
    return events


@producer(
    uses={"hcand_*"
          } | {"GenPart.*", "GenVtx.*"},
    produces={"gen_lep.*"},
    exposed=False,
)
def gen_lep_fields(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    ch_name = self.config_inst.channels.names()[0]
    hcand = events[f'hcand_{ch_name}']
    gentau = {}
    # Get index of the gen muon/electron

    # Get the the 3-vector of gen vertex
    gen_pv = get_vec_p3(events.GenVtx)
    for the_lep in ['lep0', 'lep1']:
        # Get index of the gen object that corresponds to the object from the selected dilepton pair
        gen_lep_idx = hcand[the_lep].genPartIdx
        # Get the gen object itself
        gen_lep = events.GenPart[gen_lep_idx]
        # Get secondary vertex of tau->lep decay
        gen_sv = get_vec_p3(gen_lep, 'v')
        gen_lep_vec = get_lep_p4(gen_lep).to_3D().to_pxpypz()
        # Coffea unit method does not work, so I implemented my own method
        def unit(v): return ak.where(v.mag > 0, v/v.mag, v/1.)
        gen_lep_dir = unit(gen_lep_vec)
        # Estimate tau trajectory by substracting coordinates of primary vertex from secondary vertex
        tau_trajectory = gen_pv - gen_sv
        # Calculate component of tau 3-vector, that is parallel to the decay product direction
        tau_proj = tau_trajectory.dot(gen_lep_dir)
        # Calculate Impact parameter vector as the normal component of the tau trajectory to the decay product flight direction
        gen_ip_vec = tau_trajectory - tau_proj * gen_lep_dir
        # Check for nans
        nan_mask = np.isnan(gen_ip_vec.mag)
        if ak.any(nan_mask):
            logger.warning(
                f"Obtained nan values while calculating gen impact parameter variable! Imputing them with EMPTY_FLOAT"
            )
            gen_ip_vec = ak.where(nan_mask,
                                  EMPTY_FLOAT*ak.ones_like(gen_ip_vec),
                                  gen_ip_vec)
        gentau[the_lep] = gen_lep
        for the_comp in ['x', 'y', 'z']:
            gentau[the_lep][f'gen_ip_{the_comp}'] = gen_ip_vec[the_comp]
    events = set_ak_column(events, f"gen_lep",  ak.zip(gentau))
    return events
