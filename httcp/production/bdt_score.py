# coding: utf-8

"""
Producers for the TauTheDifference BDT for signal vs. background separation of the Higgs CP analysis.
"""
import law
import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, dev_sandbox
from law.util import DotDict, InsertableDict
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
        "hcand*", "PuppiMET*", "lead_jet*", "sublead_jet*", "dijet*",
        "n_jets", "N_b_jets",
    },
    produces={
        f"bdt_raw_score_{the_output}"
        for the_output in ["gtau", "higgs", "fake"]
    } | {"bdt_cat"},
    sandbox=dev_sandbox("bash::$HTTCP_BASE/sandboxes/venv_columnar_xgb.sh"),
)
def hcp_bdt_score(
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
    features_dict = {
        'pt_1': flat_np_view(hcand.lep0.pt, axis=1),
        'pt_2': flat_np_view(hcand.lep1.pt, axis=1),
        'abs_eta_1': abs(flat_np_view(hcand.lep0.eta, axis=1)),
        'met_pt': flat_np_view(met.pt),
        'met_dphi_1': flat_np_view(hcand.delta_phi_0_met, axis=1),
        'met_dphi_2': flat_np_view(hcand.delta_phi_1_met, axis=1),
        'dR': flat_np_view(hcand.delta_r),
        'dphi': flat_np_view(hcand.delta_phi),
        'pt_tt': flat_np_view(hcand.pt_tt),
        'm_vis': flat_np_view(hcand.mass),
        'pt_vis': flat_np_view(hcand.pt_vis),
        'FastMTT_mass': flat_np_view(hcand.fastMTT.mass),
        'mt_1': flat_np_view(hcand.mt_0),
        'mt_2': flat_np_view(hcand.mt_1),
        'mt_lep': flat_np_view(hcand.mt_vis),
        'mt_tot': flat_np_view(hcand.mt_tot),
        'jpt_1': flat_np_view(events.lead_jet.pt),
        'jpt_2': flat_np_view(events.sublead_jet.pt),
        'jeta_1': flat_np_view(events.lead_jet.eta),
        'jeta_2': flat_np_view(events.sublead_jet.eta),
        'mjj': flat_np_view(events.dijet.mass),
        'jdeta': flat_np_view(events.dijet.deltaeta),
        'dijetpt': flat_np_view(events.dijet.pt),
        'n_jets': flat_np_view(events.n_jets),
        #'n_bjets': flat_np_view(events.N_b_jets),
    }
    features = pd.DataFrame.from_dict(features_dict)
    features.index = events.event
    features_even = features[features.index % 2 == 0]
    features_odd = features[features.index % 2 == 1]

    event_n = flat_np_view(events.event)

    output = np.zeros((len(event_n), 3), dtype=np.float32)
    mask_even = (event_n % 2 == 0)
    with self.evaluator:
        res_even = np.array(self.evaluator("hcp_even", features_even))
        res_odd = np.array(self.evaluator("hcp_odd", features_odd))
        output[mask_even, :] = res_even[:, :]
        output[~mask_even, :] = res_odd[:, :]
        bdt_cat = np.argmax(output, axis=1)
        for idx, the_col in enumerate(['gtau', 'higgs', 'fake']):
            events = set_ak_column_f32(
                events, f"bdt_raw_score_{the_col}", np.ascontiguousarray(output[:, idx]))
        events = set_ak_column_i32(
            events, "bdt_cat", np.ascontiguousarray(bdt_cat))

    return events


@hcp_bdt_score.requires
def hcp_bdt_score_requires(self: Producer, task: law.Task, reqs: dict) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)


@hcp_bdt_score.setup
def hcp_bdt_score_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
    task: law.task
) -> None:
    """
    Sets up XGBoost model for signal vs. bkg classification in higgs cp analysis
    """
    from httcp.ml.xgb_evaluator import XGBEvaluator

    # unpack the external files bundle and setup the evaluator
    bundle = reqs["external_files"]
    self.evaluator = XGBEvaluator()
    self.evaluator.add_model("hcp_even", bundle.files.ml_model_even.path)
    self.evaluator.add_model("hcp_odd", bundle.files.ml_model_odd.path)


@hcp_bdt_score.teardown
def hcp_bdt_score_teardown(self: Producer, **kwargs) -> None:
    """
    Stops the XGB evaluator.
    """
    if (evaluator := getattr(self, "evaluator", None)) is not None:
        evaluator.stop()
