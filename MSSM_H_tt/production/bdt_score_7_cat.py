# coding: utf-8

"""
Producers for the TauTheDifference BDT for signal vs. background separation of the Higgs MSSM analysis.
Now supports multiple mass points: one even/odd model per mass, producing per-mass outputs.

Updated for 7-class training:
  0: ggH_tautau  -> "ggh"
  1: bbH_tautau  -> "bbh"
  2: DY          -> "dy"
  3: TT          -> "tt"
  4: W+Jets      -> "wj"
  5: SingleTop   -> "st"
  6: MultiBoson  -> "mb"
"""
import law
import luigi
import functools
from columnflow.util import maybe_import, DotDict, dev_sandbox
from law.util import InsertableDict
from columnflow.columnar_util import (
    sorted_indices_from_mask,
    set_ak_column,
    has_ak_column,
    EMPTY_FLOAT,
    Route,
    flat_np_view,
    optional_column as optional,
)
from columnflow.types import Any
from columnflow.production import Producer, producer
logger = law.logger.get_logger(__name__)

np = maybe_import("numpy")
ak = maybe_import("awkward")
pd = maybe_import("pandas")
warn = maybe_import("warnings")
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

# -------------------------------------------------------------------------
# Mass points and model locations
# -------------------------------------------------------------------------
from MSSM_H_tt.config.mass_points import read_bdt_masses
MASS_POINTS = read_bdt_masses()

BDT_EOS_BASE = "/eos/user/j/jmalvaso/MSSM_H_tt_store/bdt_runs/"

def _even_path(mass: int) -> str:
    return f"{BDT_EOS_BASE}M{mass}/bst_model_M{mass}_even.json"

def _odd_path(mass: int) -> str:
    return f"{BDT_EOS_BASE}M{mass}/bst_model_M{mass}_odd.json"

# class labels in the order of training labels y = 0..6
BDT_LABELS = [
    "ggh",  # 0: ggH_tautau
    "bbh",  # 1: bbH_tautau
    "dy",   # 2: DY
    "tt",   # 3: TT
    "wj",   # 4: W+Jets
    "st",   # 5: SingleTop
    "mb",   # 6: MultiBoson
]

# list of feature names and order used in training (excluding "event")
BDT_FEATURES = [
    "D_zeta",
    "pt_1",
    "pt_2",
    "abs_eta_1",
    "abs_eta_2",
    "met_pt",
    "dR",
    "m_vis",
    "mt_tot",
    "mt_e",
    "mt_mu",
    "mt_emu",
    "ip_sig_e",
    "ip_sig_mu",
    "jpt_1",
    "jpt_2",
    "jeta_1",
    "jeta_2",
    "n_jets",
    "n_bjets",
    "pt_H",
    "fastMTT",
]

# -------------------------------------------------------------------------
# Producer
# -------------------------------------------------------------------------
@producer(
    uses={
        "event",
        "hcand*", "PuppiMET*", "lead_jet.*", "sublead_jet.*",
        "n_jets", "N_b_jets",
    },
    produces=(
        {
            f"bdt_raw_score_{lbl}_M{m}"
            for m in MASS_POINTS
            for lbl in BDT_LABELS
        }
        | {f"bdt_cat_M{m}" for m in MASS_POINTS}
    ),
    sandbox=dev_sandbox("bash::$HTTCP_BASE/sandboxes/venv_columnar_xgb.sh"),
)
def mssm_bdt_score(
    self: Producer,
    events: ak.Array,
    **kwargs,
) -> ak.Array:
    """
    Returns per-mass 7-class BDT scores and per-mass winning-category indices for each Higgs candidate.
    """
    ch_str = self.config_inst.channels.names()[0]

    hcand = events[f"hcand_{ch_str}"]
    met = events.PuppiMET

    # Event numbers (used for parity split)
    event_n = flat_np_view(events.event)

    # Lepton and MET components for features
    pt_mu   = flat_np_view(hcand.lep0.pt, axis=1)
    phi_mu  = flat_np_view(hcand.lep0.phi, axis=1)
    pt_ele  = flat_np_view(hcand.lep1.pt, axis=1)
    phi_ele = flat_np_view(hcand.lep1.phi, axis=1)
    pt_met  = flat_np_view(met.pt)
    phi_met = flat_np_view(met.phi)

    # Build features in the same way/order as in training (see get_features in training script)
    features_dict = {
        "D_zeta":     flat_np_view(events.D_zeta),
        "pt_1":       pt_mu,
        "pt_2":       pt_ele,
        "abs_eta_1":  np.abs(flat_np_view(hcand.lep0.eta, axis=1)),
        "abs_eta_2":  np.abs(flat_np_view(hcand.lep1.eta, axis=1)),
        "met_pt":     pt_met,
        "dR":         flat_np_view(hcand.delta_r),
        "m_vis":      flat_np_view(hcand.mass),
        "mt_tot":     flat_np_view(hcand.mt_tot),
        "mt_e":       flat_np_view(hcand.mt_e),
        "mt_mu":      flat_np_view(hcand.mt_mu),
        "mt_emu":     flat_np_view(hcand.mt_emu),
        "ip_sig_e":   flat_np_view(hcand.lep0.ip_sig, axis=1),
        "ip_sig_mu":  flat_np_view(hcand.lep1.ip_sig, axis=1),
        "jpt_1":      flat_np_view(events.lead_jet.pt),
        "jpt_2":      flat_np_view(events.sublead_jet.pt),
        "jeta_1":     flat_np_view(events.lead_jet.eta),
        "jeta_2":     flat_np_view(events.sublead_jet.eta),
        "n_jets":     flat_np_view(events.n_jets),
        "n_bjets":    flat_np_view(events.N_b_jets),
        "pt_H":       flat_np_view(events.pt_H),
        # additional feature used in the 7-class training
        "fastMTT":    flat_np_view(hcand.fastMTT.mass, axis=1),
    }

    # DataFrame with fixed column order matching training
    features = pd.DataFrame.from_dict(features_dict)
    features = features[BDT_FEATURES]
    features.index = event_n

    # Parity split as in training: even / odd by event number
    mask_even = (event_n % 2 == 0)
    features_even = features[mask_even]
    features_odd  = features[~mask_even]

    n_classes = len(BDT_LABELS)

    # Evaluate all masses
    from MSSM_H_tt.config.mass_points import read_bdt_masses
    MASS_POINTS = read_bdt_masses()
    with self.evaluator:
        for mass in MASS_POINTS:
            # Model keys as registered in setup
            key_even = f"bdt_even_M{mass}"
            key_odd  = f"bdt_odd_M{mass}"
  
            # run the two models (even / odd)
            res_even = np.array(self.evaluator(key_even, features_even))
            res_odd  = np.array(self.evaluator(key_odd,  features_odd))

            # stitch back together in event order
            output = np.zeros((len(event_n), n_classes), dtype=np.float32)
            output[mask_even, :]  = res_even[:, :n_classes]
            output[~mask_even, :] = res_odd[:, :n_classes]

            # per-mass winning class (0..6)
            bdt_cat = np.argmax(output, axis=1).astype(np.int32)

            # write columns for this mass, one per class
            for idx, the_col in enumerate(BDT_LABELS):
                events = set_ak_column_f32(
                    events,
                    f"bdt_raw_score_{the_col}_M{mass}",
                    np.ascontiguousarray(output[:, idx]),
                )

            events = set_ak_column_i32(
                events,
                f"bdt_cat_M{mass}",
                np.ascontiguousarray(bdt_cat),
            )
    
    return events


@mssm_bdt_score.requires
def mssm_bdt_score_requires(
    self: Producer,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    **kwargs,
    ) -> None:
    """
    No external bundle needed; models are read directly from EOS in setup.
    """
    return


@mssm_bdt_score.setup
def mssm_bdt_score_setup(
    self: Producer,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    inputs: dict[str, Any],
    reader_targets: law.util.InsertableDict,
    **kwargs,
    ) -> None:
    """
    Loads one XGBoost model pair (even/odd) for each mass point.
    """
    from MSSM_H_tt.ml.xgb_evaluator import XGBEvaluator

    self.evaluator = XGBEvaluator()
    # Register all models with unique names
    for mass in MASS_POINTS:
        self.evaluator.add_model(f"bdt_even_M{mass}", _even_path(mass))
        self.evaluator.add_model(f"bdt_odd_M{mass}",  _odd_path(mass))


@mssm_bdt_score.teardown
def mssm_bdt_score_teardown(self: Producer, **kwargs) -> None:
    """
    Stops the XGB evaluator.
    """
    if (evaluator := getattr(self, "evaluator", None)) is not None:
        evaluator.stop()
