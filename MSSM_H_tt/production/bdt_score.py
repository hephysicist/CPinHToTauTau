# coding: utf-8

"""
Producers for the TauTheDifference BDT for signal vs. background separation of the Higgs MSSM analysis.
Now supports multiple mass points: one even/odd model per mass, producing per-mass outputs.
"""
import law
import luigi
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, dev_sandbox, DotDict, InsertableDict
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column, flat_np_view

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

BDT_EOS_BASE = "/eos/user/j/jmalvaso/SWAN_projects/XGBoost_MSSM/"
def _even_path(mass: int) -> str:
    return f"{BDT_EOS_BASE}M{mass}/bst_model_M{mass}_even.json"
def _odd_path(mass: int) -> str:
    return f"{BDT_EOS_BASE}M{mass}/bst_model_M{mass}_odd.json"

# -------------------------------------------------------------------------
# Producer
# -------------------------------------------------------------------------
@producer(
    uses={
        "event",
        "hcand*", "PuppiMET*", "leading_jet*", "subleading_jet*",
        "n_jets", "N_b_jets",
    },
    produces=(
        {
            f"bdt_raw_score_{lbl}_M{m}"
            for m in MASS_POINTS
            for lbl in ["sig", "dy", "tt", "wj"]
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
    Returns per-mass BDT scores and per-mass winning-category indices for each Higgs candidate.
    """
    ch_str = self.config_inst.channels.names()[0]

    hcand = events[f"hcand_{ch_str}"]
    met = events.PuppiMET

    # Lepton and MET components for features
    pt_mu   = flat_np_view(hcand.lep0.pt, axis=1)
    phi_mu  = flat_np_view(hcand.lep0.phi, axis=1)
    pt_ele  = flat_np_view(hcand.lep1.pt, axis=1)
    phi_ele = flat_np_view(hcand.lep1.phi, axis=1)
    pt_met  = flat_np_view(met.pt)
    phi_met = flat_np_view(met.phi)

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
        "ip_sig_e":   flat_np_view(hcand.lep0.ip_sig),
        "ip_sig_mu":  flat_np_view(hcand.lep1.ip_sig),
        "jpt_1":      flat_np_view(events.leading_jet_pt),
        "jpt_2":      flat_np_view(events.subleading_jet_pt),
        "jeta_1":     flat_np_view(events.leading_jet_eta),
        "jeta_2":     flat_np_view(events.subleading_jet_eta),
        "n_jets":     flat_np_view(events.n_jets),
        "n_bjets":    flat_np_view(events.N_b_jets),
        "pt_H":       flat_np_view(events.pt_H),
    }
    features = pd.DataFrame.from_dict(features_dict)
    features.index = events.event

    features_even = features[features.index % 2 == 0]
    features_odd  = features[features.index % 2 == 1]

    event_n = flat_np_view(events.event)
    mask_even = (event_n % 2 == 0)

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
            output = np.zeros((len(event_n), 4), dtype=np.float32)
            output[mask_even, :]  = res_even[:, :]
            output[~mask_even, :] = res_odd[:, :]

            # per-mass winning class (0..3)
            bdt_cat = np.argmax(output, axis=1)

            # write columns for this mass
            for idx, the_col in enumerate(["sig", "dy", "tt", "wj"]):
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
def mssm_bdt_score_requires(self: Producer, reqs: dict) -> None:
    """
    No external bundle needed; models are read directly from EOS in setup.
    """
    return


@mssm_bdt_score.setup
def mssm_bdt_score_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
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
