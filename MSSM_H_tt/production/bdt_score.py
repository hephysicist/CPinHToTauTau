# coding: utf-8

"""
Producers for the TauTheDifference BDT for signal vs. background separation of the Higgs MSSM analysis.
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


@producer(
    uses={
        "event",
        "hcand*", "PuppiMET*", "leading_jet*", "subleading_jet*",
        "n_jets", "N_b_jets",
    },
    produces={
        f"bdt_raw_score_{the_output}"
        for the_output in ["sig", "dy", "tt", "wj"]
    } | {"bdt_cat"},
    sandbox=dev_sandbox("bash::$HTTCP_BASE/sandboxes/venv_columnar_xgb.sh"),
)
def mssm_bdt_score(
    self: Producer,
    events: ak.Array,
    **kwargs,
) -> ak.Array:
    """
    Returns the bdt score per higgs candidate object. This score will be used to create signal categories.
    """

    ch_str = self.config_inst.channels.names()[0]

    hcand = events[f'hcand_{ch_str}']
    met = events.PuppiMET
    # Lepton and MET components for pt_H
    pt_mu   = flat_np_view(hcand.lep0.pt, axis=1)
    phi_mu  = flat_np_view(hcand.lep0.phi, axis=1)
    pt_ele  = flat_np_view(hcand.lep1.pt, axis=1)
    phi_ele = flat_np_view(hcand.lep1.phi, axis=1)
    pt_met  = flat_np_view(met.pt)
    phi_met = flat_np_view(met.phi)

    # Vector sum in px, py → pt_H
    px_H = pt_mu * np.cos(phi_mu) + pt_ele * np.cos(phi_ele) + pt_met * np.cos(phi_met)
    py_H = pt_mu * np.sin(phi_mu) + pt_ele * np.sin(phi_ele) + pt_met * np.sin(phi_met)
    pt_H = np.sqrt(px_H**2 + py_H**2)

    features_dict = {
        'D_zeta':     flat_np_view(events.D_zeta),
        'pt_1':       pt_mu,
        'pt_2':       pt_ele,
        'abs_eta_1':  np.abs(flat_np_view(hcand.lep0.eta, axis=1)),
        'abs_eta_2':  np.abs(flat_np_view(hcand.lep1.eta, axis=1)),
        'met_pt':     pt_met,
        'dR':         flat_np_view(hcand.delta_r),
        'm_vis':      flat_np_view(hcand.mass),
        'mt_tot':     flat_np_view(hcand.mt_tot),
        'mt_e':       flat_np_view(hcand.mt_e),
        'mt_mu':      flat_np_view(hcand.mt_mu),
        'mt_emu':     flat_np_view(hcand.mt_emu),
        'ip_sig_e':   flat_np_view(hcand.lep0.ip_sig),
        'ip_sig_mu':  flat_np_view(hcand.lep1.ip_sig),
        'jpt_1':      flat_np_view(events.leading_jet_pt),
        'jpt_2':      flat_np_view(events.subleading_jet_pt),
        'jeta_1':     flat_np_view(events.leading_jet_eta),
        'jeta_2':     flat_np_view(events.subleading_jet_eta),
        'n_jets':     flat_np_view(events.n_jets),
        'n_bjets':    flat_np_view(events.N_b_jets),
        'pt_H':       pt_H
    }
    features = pd.DataFrame.from_dict(features_dict)
    features.index = events.event
    features_even = features[features.index % 2 == 0]
    features_odd = features[features.index % 2 == 1]

    event_n = flat_np_view(events.event)

    output = np.zeros((len(event_n), 4), dtype=np.float32)
    mask_even = (event_n % 2 == 0)
    with self.evaluator:
        res_even = np.array(self.evaluator("bdt_even", features_even))
        res_odd = np.array(self.evaluator("bdt_odd", features_odd))
        output[mask_even, :] = res_even[:, :]
        output[~mask_even, :] = res_odd[:, :]
        bdt_cat = np.argmax(output, axis=1)

        for idx, the_col in enumerate(["sig", "dy", "tt", "wj"]):
            events = set_ak_column_f32(
                events, f"bdt_raw_score_{the_col}", np.ascontiguousarray(output[:, idx]))
        events = set_ak_column_i32(
            events, "bdt_cat", np.ascontiguousarray(bdt_cat))

    return events


@mssm_bdt_score.requires
def mssm_bdt_score_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)


@mssm_bdt_score.setup
def mssm_bdt_score_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    """
    Sets up XGBoost model for signal vs. bkg classification in MSSM analysis
    """
    from MSSM_H_tt.ml.xgb_evaluator import XGBEvaluator

    # unpack the external files bundle and setup the evaluator
    bundle = reqs["external_files"]
    self.evaluator = XGBEvaluator()
    mass =  self.config_inst.x.bdt_mass["mass"]
    self.evaluator.add_model("bdt_even", bundle.files[f"ml_model_even_{mass}"].path)
    self.evaluator.add_model("bdt_odd", bundle.files[f"ml_model_odd_{mass}"].path)


@mssm_bdt_score.teardown
def mssm_bdt_score_teardown(self: Producer, **kwargs) -> None:
    """
    Stops the XGB evaluator.
    """
    if (evaluator := getattr(self, "evaluator", None)) is not None:
        evaluator.stop()