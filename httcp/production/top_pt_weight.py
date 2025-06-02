# coding: utf-8

"""
Column producers related to top quark pt reweighting.
"""

import law

from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

ak = maybe_import("awkward")
np = maybe_import("numpy")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

logger = law.logger.get_logger(__name__)


@producer(
    uses={"GenPart.{pdgId,statusFlags}"},
    # requested GenPartonTop columns, passed to the *uses* and *produces*
    produced_top_columns={"pt"},
    mc_only=True,
)
def gen_parton_top(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Produce parton-level top quarks (before showering and detector simulation).
    Creates new collection named "GenPartonTop"

    *produced_top_columns* can be adapted to change the columns that will be produced
    for the GenPartonTop collection.

    The function is skipped when the dataset is data or when it does not have the tag *has_top*.

    :param events: awkward array containing events to process
    """
    # find parton-level top quarks
    abs_id = abs(events.GenPart.pdgId)
    t = events.GenPart[abs_id == 6]
    bit = 13 # It is crucial that the pT value used when you calculate the weights is derived from the 'isLastCopy' definition of the parton-level top quark (after radiation and before decay) in Pythia8
    bit_mask = 1<<(bit-1)
    Filter_Bits = ((t.statusFlags & bit_mask) != 0)
    t = t[Filter_Bits]
    t = t[~ak.is_none(t, axis=1)]

    # save the column
    events = set_ak_column(events, "GenPartonTop", t)

    return events


@gen_parton_top.init
def gen_parton_top_init(self: Producer) -> bool:
    for col in self.produced_top_columns:
        self.uses.add(f"GenPart.{col}")
        self.produces.add(f"GenPartonTop.{col}")
        

@producer(
    uses={
        "GenPartonTop.pt",
    },
    produces={
        "top_pt_weight", "top_pt_weight_up", "top_pt_weight_down",
    },
    get_top_pt_config=(lambda self: self.config_inst.x.top_pt_reweighting_params),
)
def top_pt_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Compute SF to be used for top pt reweighting.

    The *GenPartonTop.pt* column can be produced with the :py:class:`gen_parton_top` Producer.

    The SF should *only be applied in ttbar MC* as an event weight and is computed
    based on the gen-level top quark transverse momenta.

    The function is skipped when the dataset is data or when it does not have the tag *is_ttbar*.

    The top pt reweighting parameters should be given as an auxiliary entry in the config:

    .. code-block:: python

        cfg.x.top_pt_reweighting_params = {
            "a": 0.0615,
            "a_up": 0.0615 * 1.5,
            "a_down": 0.0615 * 0.5,
            "b": -0.0005,
            "b_up": -0.0005 * 1.5,
            "b_down": -0.0005 * 0.5,
        }

    *get_top_pt_config* can be adapted in a subclass in case it is stored differently in the config.

    :param events: awkward array containing events to process
    """

    # get SF function parameters from config
    params = self.get_top_pt_config()

    # check the number of gen tops
    if ak.any(ak.num(events.GenPartonTop, axis=1) != 2):
        logger.warning("There are events with != 2 GenPartonTops. This producer should only run for ttbar")

    # clamp top pT < 500 GeV
    pt_clamped = ak.where(events.GenPartonTop.pt > 500.0, 500.0, events.GenPartonTop.pt)
    for variation in ("", "_up", "_down"):
        # evaluate SF function
        sf = np.exp(params[f"a{variation}"] + params[f"b{variation}"] * pt_clamped)
        extrapolation_sf_13p6 = 0.991 + 0.000075 * pt_clamped #Extrapolation suggested by top group to go from 13 TeV to 13.6 TeV
        sf = sf*extrapolation_sf_13p6
        # compute weight from SF product for top and anti-top
        weight = np.sqrt(np.prod(sf, axis=1))

        # write out weights
        events = set_ak_column(events, f"top_pt_weight{variation}", ak.fill_none(weight, 1.0))

    return events

