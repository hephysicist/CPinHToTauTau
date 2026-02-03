"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT
from columnflow.util import maybe_import, DotDict
from MSSM_H_tt.util import get_lep_p4, get_vec_p3, to_pt_eta_phi_m

np = maybe_import("numpy")
ak = maybe_import("awkward")

import functools
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


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
    #Create a template array filled with zeros
    channel_id  = ak.zeros_like(ak.local_index(events.event), dtype = np.uint8)
    
    #Create a mask array to check for channel orthogonality 
    and_mask = ak.ones_like(ak.local_index(events.event), dtype = np.bool_)
    
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
        "jet_pt_30": sorted_jets.pt > 30.0,
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
          ["pt", "eta", "phi", "mass", "btagDeepFlavB", "pass_tightID_lep_veto"]} | {"hcand_*"},
    produces={
        "N_b_jets",
        "lead_b_jet.*",
        "sublead_b_jet.*",
        "di_b_jet.*",
    },
    exposed=False,
)
def number_b_jet(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    Produces:
      - N_b_jets
      - lead_b_jet.{pt,eta,phi,mass}
      - sublead_b_jet.{pt,eta,phi,mass} (and alias sublad_b_jet.*)
      - di_b_jet.{pt,eta,phi,mass,deltaeta,deltaphi,delta_r}
    """
    year = self.config_inst.x.year
    tag = self.config_inst.x.tag
    btag_wp = self.config_inst.x.btag_working_points[year][tag].deepjet.medium

    # sort jets by pt
    jet_pt_sorted_idx = ak.argsort(events.Jet.pt, axis=1, ascending=False)
    sorted_jets = events.Jet[jet_pt_sorted_idx]

    # base (b-jet) selection
    jet_selections = {
        "jet_pt_20": sorted_jets.pt > 20.0,
        "jet_eta_2.5": abs(sorted_jets.eta) < 2.5,
        "jet_id": sorted_jets.pass_tightID_lep_veto,
        "btag_wp_medium": sorted_jets.btagDeepFlavB >= btag_wp,
    }
    jet_obj_mask = ak.ones_like(jet_pt_sorted_idx, dtype=np.bool_)
    for the_sel in jet_selections.values():
        jet_obj_mask = jet_obj_mask & the_sel

    # use first configured channel (same pattern as jet_pt_def)
    ch_str = self.config_inst.channels.names()[0]
    hcand = events[f"hcand_{ch_str}"]

    # Δφ wrapping helper (vectorized)
    def _wrap_delta_phi(dphi):
        dphi = ak.where(dphi > np.pi, dphi - 2.0 * np.pi, dphi)
        dphi = ak.where(dphi < -np.pi, dphi + 2.0 * np.pi, dphi)
        return dphi

    # remove jets matched to leptons (seed jetIdx) and require ΔR(jet,lep) > 0.4 for both leptons
    mask = jet_obj_mask
    for lep_str in ["lep0", "lep1"]:
        lep = ak.firsts(hcand[lep_str])
        seed_idx = ak.fill_none(lep.jetIdx, -1)

        # veto the seed jet (same indexing space as sorted jets)
        mask = mask & (jet_pt_sorted_idx != seed_idx)

        # ΔR cleaning
        dphi = _wrap_delta_phi(sorted_jets.phi - lep.phi)
        deta = sorted_jets.eta - lep.eta
        dr = np.sqrt(dphi**2 + deta**2)
        mask = mask & ak.fill_none((dr > 0.4), False)

    # selected b-jets after cleaning
    sel_bjets = ak.drop_none(ak.mask(sorted_jets, mask))
    nbjets = ak.num(sel_bjets, axis=1)

    # build leading/subleading and dijet (using same utilities as jet_pt_def)
    empty_p4 = ak.zeros_like(get_lep_p4(hcand.lep0))
    sel_bjets_p4 = get_lep_p4(sel_bjets)

    lead_bjet_p4 = ak.where(nbjets > 0, sel_bjets_p4[:, :1], empty_p4)
    sublead_bjet_p4 = ak.where(nbjets > 1, sel_bjets_p4[:, 1:2], empty_p4)

    di_bjet_p4 = to_pt_eta_phi_m(lead_bjet_p4 + sublead_bjet_p4)

    # fill output records with EMPTY_FLOAT when not present
    empty_float = lambda arr: ak.full_like(arr, EMPTY_FLOAT)

    lead_b_jet = {}
    lead_mask = (lead_bjet_p4.pt > 20)
    sublead_b_jet = {}
    sublead_mask = (sublead_bjet_p4.pt > 20)
    di_b_jet = {}
    di_mask = (di_bjet_p4.pt > 20)

    for var in ["pt", "eta", "phi", "mass"]:
        lead_b_jet[var] = ak.where(
            lead_mask, getattr(lead_bjet_p4, var), empty_float(getattr(lead_bjet_p4, var))
        )
        sublead_b_jet[var] = ak.where(
            sublead_mask, getattr(sublead_bjet_p4, var), empty_float(getattr(sublead_bjet_p4, var))
        )
        di_b_jet[var] = ak.where(
            di_mask, getattr(di_bjet_p4, var), empty_float(getattr(di_bjet_p4, var))
        )

    # dijet angular separations (only meaningful if >=2 b-jets)
    for func_name in ["deltaeta", "deltaphi", "delta_r"]:
        func = getattr(lead_bjet_p4, func_name)
        di_b_jet[func_name] = ak.where(
            nbjets >= 2, func(sublead_bjet_p4), empty_float(di_bjet_p4.pt)
        )

    # write columns
    events = set_ak_column(events, "lead_b_jet", ak.zip(lead_b_jet))
    events = set_ak_column(events, "sublead_b_jet", ak.zip(sublead_b_jet))
    events = set_ak_column(events, "di_b_jet", ak.zip(di_b_jet))
    events = set_ak_column(events, "N_b_jets", nbjets)
    return events