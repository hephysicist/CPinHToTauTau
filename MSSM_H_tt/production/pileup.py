# coding: utf-8

"""
Column production methods related to pileup weights.
"""

import functools

import law

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, DotDict
from columnflow.columnar_util import set_ak_column
from columnflow.types import Any

np = maybe_import("numpy")
ak = maybe_import("awkward")


# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

logger = law.logger.get_logger(__name__)


@producer(
    uses={"Pileup.nTrueInt"},
    produces={"pu_weight{,_up,_down}"},
    # only run on mc
    mc_only=True,
    # function to determine the correction file
    get_pileup_file=(lambda self, external_files: external_files.pu_sf),
)
def pu_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Based on the number of primary vertices, assigns each event pileup weights using correctionlib.
    """
    # map the variable names from the corrector to our columns
    variable_map = {
        "NumTrueInteractions": events.Pileup.nTrueInt,
    }

    rename_systs = {"pu_weight" :'pu_weight',
                    "pu_weight_minbias_xs_up"  : 'pu_weight_up',
                    "pu_weight_minbias_xs_down": 'pu_weight_down'
    }   
    for column_name, syst in (
        ("pu_weight", "nominal"),
        ("pu_weight_minbias_xs_up", "up"),
        ("pu_weight_minbias_xs_down", "down"),
    ):
        # get the inputs for this type of variation
        variable_map_syst = {**variable_map, "weights": syst}
        inputs = [variable_map_syst[inp.name] for inp in self.pileup_corrector.inputs]
        # evaluate and store the produced column
        pu_weight = self.pileup_corrector.evaluate(*inputs)
        pu_weight[pu_weight > 300] = 0 # remove events with abnormally high pileup weights (should be few per dataset)
        events = set_ak_column(events, rename_systs[column_name], pu_weight, value_type=np.float32)
    
    return events


@pu_weight.requires
def pu_weight_requires(
    self: Producer,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    **kwargs,
) -> None:
    """
    Adds the requirements needed the underlying task to derive the pileup weights into *reqs*.
    """
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)


@pu_weight.setup
def pu_weight_setup(
    self: Producer,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    inputs: dict[str, Any],
    reader_targets: law.util.InsertableDict,
    **kwargs,
) -> None:
    """
    Loads the pileup calculator from the external files bundle and saves them in the
    py:attr:`pileup_corrector` attribute for simpler access in the actual callable.
    """
    bundle = reqs["external_files"]

    # create the corrector
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_pileup_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )

    # check
    if len(correction_set.keys()) != 1:
        raise Exception("Expected exactly one type of pileup correction")

    corrector_name = list(correction_set.keys())[0]
    self.pileup_corrector = correction_set[corrector_name]
