# coding: utf-8

"""
Jet energy corrections and jet resolution smearing.
"""
from pprint import pprint

import functools

from columnflow.types import Any
from columnflow.calibration import Calibrator, calibrator
from httcp.calibration.util import ak_random, propagate_met
from columnflow.production.util import attach_coffea_behavior
from columnflow.util import maybe_import
from law.util import InsertableDict, DotDict
from columnflow.columnar_util import set_ak_column, layout_ak_array, optional_column as optional
import law

np = maybe_import("numpy")
ak = maybe_import("awkward")
correctionlib = maybe_import("correctionlib")


#
# helper functions
#

set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


import difflib

def get_evaluators(
    correction_set: correctionlib.highlevel.CorrectionSet,
    names: list[str],
) -> list[Any]:
    """
    Helper function to get a list of correction evaluators from a
    :external+correctionlib:py:class:`correctionlib.highlevel.CorrectionSet` object given
    a list of *names*. The *names* can refer to either simple or compound
    corrections.

    :param correction_set: evaluator provided by :external+correctionlib:doc:`index`
    :param names: List of names of corrections to be applied
    :raises RuntimeError: If a requested correction in *names* is not available
    :return: List of compounded corrections, see
        :external+correctionlib:py:class:`correctionlib.highlevel.CorrectionSet`
    """
    available_keys = set(correction_set.keys()).union(correction_set.compound.keys())
    corrected_names = []
    for name in names:
        if name not in available_keys:
            # Find the closest match using difflib
            closest_matches = difflib.get_close_matches(name, available_keys, n=1)
            if closest_matches:
                closest_match = closest_matches[0]
                print(
                    f"Correction '{name}' not found. Using closest match: '{closest_match}'",
                )
                corrected_names.append(closest_match)
            else:
                raise RuntimeError(f"Correction '{name}' not found and no close match available.")
        else:
            corrected_names.append(name)
    
    # Retrieve the evaluators
    return [
        correction_set.compound[name]
        if name in correction_set.compound
        else correction_set[name]
        for name in corrected_names
    ]

def ak_evaluate(evaluator: correctionlib.highlevel.Correction, *args) -> float:
    """
    Evaluate a :external+correctionlib:py:class:`correctionlib.highlevel.Correction`
    using one or more :external+ak:py:class:`awkward arrays <ak.Array>` as inputs.

    :param evaluator: Evaluator instance
    :raises ValueError: If no :external+ak:py:class:`awkward arrays <ak.Array>` are provided
    :return: The correction factor derived from the input arrays
    """
    # fail if no arguments
    if not args:
        raise ValueError("Expected at least one argument.")

    # collect arguments that are awkward arrays
    ak_args = [
        arg for arg in args if isinstance(arg, ak.Array)
    ]

    # broadcast akward arrays together and flatten
    if ak_args:
        bc_args = ak.broadcast_arrays(*ak_args)
        flat_args = (
            np.asarray(ak.flatten(bc_arg, axis=None))
            for bc_arg in bc_args
        )
        output_layout_array = bc_args[0]
    else:
        flat_args = iter(())
        output_layout_array = None

    # multiplex flattened and non-awkward inputs
    all_flat_args = [
        next(flat_args) if isinstance(arg, ak.Array) else arg
        for arg in args
    ]

    # apply evaluator to flattened/multiplexed inputs
    result = evaluator.evaluate(*all_flat_args)

    # apply broadcasted layout to result
    if output_layout_array is not None:
        result = layout_ak_array(result, output_layout_array)

    return result


#
# jet energy corrections
#

# define default functions for jec calibrator
def get_jec_file_default(self, external_files: DotDict) -> str:
    """
    Function to obtain external jec files.

    By default, this function extracts the location of the jec correction
    files from the current config instance *config_inst*:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "jet_jerc": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/JME/2017_UL/jet_jerc.json.gz",
        })

    :param external_files: Dictionary containing the information about the file location
    :return: path or url to correction file(s)
    """ # noqa
    return external_files.jet_jerc


def get_jec_config_default(self) -> DotDict:
    """Load config relevant to the JEC corrections.

    By default, this is extracted from the current *config_inst*:

    .. code-block:: python

        self.config_inst.x.jec

    Used in :py:meth:`~.jec.setup_func`.

    :return: Dictionary containing configurations for JEC callibrations
    """
    return self.config_inst.x.jec


@calibrator(
    uses={
        optional("run"),
        "Jet.pt", "Jet.eta", "Jet.phi", "Jet.mass", "Jet.area", "Jet.rawFactor",
        "Jet.jetId",
        optional("fixedGridRhoFastjetAll"),
        optional("Rho.fixedGridRhoFastjetAll"),
        attach_coffea_behavior,
    },
    produces={
        "Jet.pt", "Jet.phi", "Jet.mass", "Jet.rawFactor",
    },
    # custom uncertainty sources, defaults to config when empty
    uncertainty_sources=None,
    # toggle for propagation to PuppiMET
    propagate_met=True,
    # # function to determine the correction file
    get_jec_file=get_jec_file_default,
    # # function to determine the jec configuration dict
    get_jec_config=get_jec_config_default,
    # function to update variables before jec corrector call
    update_corrector_variables=(lambda self, corrector, variables: variables),
)

def jec(
    self: Calibrator,
    events: ak.Array,
    min_pt_met_prop: float = 15.0,
    max_eta_met_prop: float = 5.2,
    **kwargs,
) -> ak.Array:
    """Performs the jet energy corrections and uncertainty shifts using the
    :external+correctionlib:doc:`index`, optionally
    propagating the changes to the PuppiMET.

    Requires an external file in the config under ``jet_jerc``:

    .. code-block:: python

        cfg.x.external_files = DotDict.wrap({
            "jet_jerc": "/afs/cern.ch/work/m/mrieger/public/mirrors/jsonpog-integration-9ea86c4c/POG/JME/2017_UL/jet_jerc.json.gz",
        })

    *get_jec_file* can be adapted in a subclass in case it is stored differently in the
    external files

    The jec configuration should be an auxiliary entry in the config, specifying the correction
    details under "jec":

    .. code-block:: python

        cfg.x.jec = {
            "campaign": "Summer19UL17",
            "version": "V5",
            "jet_type": "AK4PFchs",
            "levels": ["L1L2L3Res"],  # or individual correction levels
            "levels_for_type1_met": ["L1FastJet"],
            "uncertainty_sources": [
                "Total",
                "CorrelationGroupMPFInSitu",
                "CorrelationGroupIntercalibration",
                "CorrelationGroupbJES",
                "CorrelationGroupFlavor",
                "CorrelationGroupUncorrelated",
            ]
        }

    *get_jec_config* can be adapted in a subclass in case it is stored differently in the config.

    If running on data, the datasets must have an auxiliary field *jec_era* defined, e.g. "RunF",
    or an auxiliary field *era*, e.g. "F".

    This instance of :py:class:`~columnflow.calibration.Calibrator` is
    initialized with the following parameters by default:

    :param events: awkward array containing events to process

    :param min_pt_met_prop: If *propagate_met* variable is ``True`` propagate the updated jet values
        to the missing transverse energy (PuppiMET) using
        :py:func:`~columnflow.calibration.util.propagate_met` for events where
        ``met.pt > *min_pt_met_prop*``.
    :param max_eta_met_prop: If *propagate_met* variable is ``True`` propagate the updated jet
        values to the missing transverse energy (PuppiMET) using
        :py:func:`~columnflow.calibration.util.propagate_met` for events where
        ``met.eta > *min_eta_met_prop*``.
    """ # noqa
    # calculate uncorrected pt, mass
    events = set_ak_column_f32(events, "Jet.pt_raw", events.Jet.pt * (1 - events.Jet.rawFactor))
    events = set_ak_column_f32(events, "Jet.mass_raw", events.Jet.mass * (1 - events.Jet.rawFactor))
    
    def correct_jets(pt, phi, area, eta, rho,run, evaluator_key="jec"):
        # variable naming convention
        variable_map = {
            "JetA": area,
            "JetEta": eta,
            "JetPt": pt,
            "JetPhi": phi,
            "Rho": ak.values_astype(rho, np.float32),
            "run": run,
        }

        # apply all correctors sequentially, updating the pt each time
        full_correction = ak.ones_like(pt, dtype=np.float32)
        

        for corrector in self.evaluators[evaluator_key]:
            # determine correct inputs (change depending on corrector)
            inputs = [
                variable_map[inp.name]
                for inp in corrector.inputs
            ]
            correction = ak_evaluate(corrector, *inputs)
            # update pt for subsequent correctors
            #pprint(corrector.__dict__)  # If `corrector` is a custom object with attributes
            variable_map["JetPt"] = variable_map["JetPt"] * correction
            full_correction = full_correction * correction

        return full_correction

    # obtain rho, which might be located at different routes, depending on the nano version
    rho = (
        events.fixedGridRhoFastjetAll
        if "fixedGridRhoFastjetAll" in events.fields else
        events.Rho.fixedGridRhoFastjetAll
    )

    # correct jets with only a subset of correction levels
    # (for calculating TypeI PuppiMET correction)
    if self.propagate_met:
        # get correction factors
        jec_factors_subset_type1_met = correct_jets(
            pt=events.Jet.pt_raw,
            phi=events.Jet.phi,
            eta=events.Jet.eta,
            area=events.Jet.area,
            rho=rho,
            run=events.run,
            evaluator_key="jec_subset_type1_met",
        )
        
        # temporarily apply the new factors with only subset of corrections
        events = set_ak_column_f32(events, "Jet.pt", events.Jet.pt_raw * jec_factors_subset_type1_met)
        events = set_ak_column_f32(events, "Jet.mass", events.Jet.mass_raw * jec_factors_subset_type1_met)
        events = self[attach_coffea_behavior](events, collections=["Jet"], **kwargs)

        # store pt and phi of the full jet system for PuppiMET propagation, including a selection in raw info
        # see https://twiki.cern.ch/twiki/bin/view/CMS/JECAnalysesRecommendations?rev=19#Minimum_jet_selection_cuts
        met_prop_mask = (events.Jet.pt_raw > min_pt_met_prop) & (abs(events.Jet.eta) < max_eta_met_prop)
        jetsum = events.Jet[met_prop_mask].sum(axis=1)
        jetsum_pt_subset_type1_met = jetsum.pt
        jetsum_phi_subset_type1_met = jetsum.phi

    # factors for full jet correction with all levels
    jec_factors = correct_jets(
        pt=events.Jet.pt_raw,
        phi=events.Jet.phi,
        eta=events.Jet.eta,
        area=events.Jet.area,
        rho=rho,
        run=events.run,
        evaluator_key="jec",
    )
    # map jet phi into [-pi, pi] 
    jet_phi_nc = events.Jet.phi
    jet_phi = ak.where(
            np.abs(jet_phi_nc) > np.pi,
            jet_phi_nc - 2 * np.pi * np.sign(jet_phi_nc),
            jet_phi_nc)
    events = set_ak_column_f32(events, "Jet.phi", jet_phi)
    
    # apply full jet correction
    events = set_ak_column_f32(events, "Jet.pt", events.Jet.pt_raw * jec_factors)
    events = set_ak_column_f32(events, "Jet.mass", events.Jet.mass_raw * jec_factors)
    rawFactor = ak.nan_to_num(1 - events.Jet.pt_raw / events.Jet.pt, nan=0.0)
    events = set_ak_column_f32(events, "Jet.rawFactor", rawFactor)
    events = self[attach_coffea_behavior](events, collections=["Jet"], **kwargs)

    # nominal met propagation
    if self.propagate_met:
        # get pt and phi of all jets after correcting
        jetsum = events.Jet[met_prop_mask].sum(axis=1)
        jetsum_pt_all_levels = jetsum.pt
        jetsum_phi_all_levels = jetsum.phi

        # propagate changes to PuppiMET, starting from jets corrected with subset of JEC levels
        # (recommendation is to propagate only L2 corrections and onwards)

        met_pt, met_phi = propagate_met(
            jetsum_pt_subset_type1_met,
            jetsum_phi_subset_type1_met,
            jetsum_pt_all_levels,
            jetsum_phi_all_levels,
            events.RawPuppiMET.pt,
            events.RawPuppiMET.phi,
        )
        events = set_ak_column_f32(events, "PuppiMET.pt", met_pt)
        events = set_ak_column_f32(events, "PuppiMET.phi", met_phi)

    # variable naming conventions
    variable_map = {
        "JetEta": events.Jet.eta,
        "JetPt": events.Jet.pt_raw,
        "JetPhi": events.Jet.phi,
    }

    # jet energy uncertainty components
    for name, evaluator in self.evaluators["junc"].items():
        # get uncertainty
        inputs = [variable_map[inp.name] for inp in evaluator.inputs]
        jec_uncertainty = ak_evaluate(evaluator, *inputs)

        # apply jet uncertainty shifts
        events = set_ak_column_f32(events, f"Jet.pt_jec_{name}_up", events.Jet.pt * (1.0 + jec_uncertainty))
        events = set_ak_column_f32(events, f"Jet.pt_jec_{name}_down", events.Jet.pt * (1.0 - jec_uncertainty))
        events = set_ak_column_f32(events, f"Jet.mass_jec_{name}_up", events.Jet.mass * (1.0 + jec_uncertainty))
        events = set_ak_column_f32(events, f"Jet.mass_jec_{name}_down", events.Jet.mass * (1.0 - jec_uncertainty))
        
        
        # propagate shifts to PuppiMET
        if self.propagate_met:
            jet_pt_up = events.Jet[met_prop_mask][f"pt_jec_{name}_up"]
            jet_pt_down = events.Jet[met_prop_mask][f"pt_jec_{name}_down"]
            met_pt_up, met_phi_up = propagate_met(
                jetsum_pt_all_levels,
                jetsum_phi_all_levels,
                jet_pt_up,
                events.Jet[met_prop_mask].phi,
                met_pt,
                met_phi,
            )
            met_pt_down, met_phi_down = propagate_met(
                jetsum_pt_all_levels,
                jetsum_phi_all_levels,
                jet_pt_down,
                events.Jet[met_prop_mask].phi,
                met_pt,
                met_phi,
            )
            
            events = set_ak_column_f32(events, f"PuppiMET.pt_jec_{name}_up", met_pt_up)
            events = set_ak_column_f32(events, f"PuppiMET.pt_jec_{name}_down", met_pt_down)
            events = set_ak_column_f32(events, f"PuppiMET.phi_jec_{name}_up", met_phi_up)
            events = set_ak_column_f32(events, f"PuppiMET.phi_jec_{name}_down", met_phi_down)

    return events


@jec.init
def jec_init(self: Calibrator) -> None:
    jec_cfg = self.get_jec_config()

    sources = self.uncertainty_sources
    if sources is None:
        sources = jec_cfg.uncertainty_sources

    # add shifted jet variables
    self.produces |= {
        f"Jet.{shifted_var}_jec_{junc_name}_{junc_dir}"
        for shifted_var in ("pt", "mass")
        for junc_name in sources
        for junc_dir in ("up", "down")
    }
    self.uses.add("run")
    # add PuppiMET variables
    if self.propagate_met:
        self.uses |= {"RawPuppiMET.pt", "RawPuppiMET.phi","PuppiMET.pt", "PuppiMET.phi"}
        self.produces |= {"PuppiMET.pt", "PuppiMET.phi"}

        # add shifted PuppiMET variables
        self.produces |= {
            f"PuppiMET.{shifted_var}_jec_{junc_name}_{junc_dir}"
            for shifted_var in ("pt", "phi")
            for junc_name in sources
            for junc_dir in ("up", "down")
        }


@jec.requires
def jec_requires(
    self: Calibrator,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    **kwargs,
    ) -> None:
    if "external_files" in reqs:
        return

    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(task)


@jec.setup
def jec_setup(
    self: Calibrator,
    task: law.Task,
    reqs: dict[str, DotDict[str, Any]],
    inputs: dict[str, Any],
    reader_targets: law.util.InsertableDict,
    **kwargs,
    ) -> None:
    """
    Load the correct jec files using the :py:func:`from_string` method of the
    :external+correctionlib:py:class:`correctionlib.highlevel.CorrectionSet`
    function and apply the corrections as needed.

    The source files for the :external+correctionlib:py:class:`correctionlib.highlevel.CorrectionSet`
    instance are extracted with the :py:meth:`~.jec.get_jec_file`.

    Uses the member function :py:meth:`~.jec.get_jec_config` to construct the
    required keys, which are based on the following information about the JEC:

        - levels
        - campaign
        - version
        - jet_type

    A corresponding example snippet wihtin the *config_inst* could like something
    like this:

    .. code-block:: python

        cfg.x.jec = DotDict.wrap({
            # campaign name for this JEC correctiono
            "campaign": f"Summer19UL{year2}{jerc_postfix}",
            # version of the corrections
            "version": "V7",
            # Type of jets that the corrections should be applied on
            "jet_type": "AK4PFchs",
            # relevant levels in the derivation process of the JEC
            "levels": ["L1FastJet", "L2Relative", "L2L3Residual", "L3Absolute"],
            # relevant levels in the derivation process of the Type 1 PuppiMET JEC
            "levels_for_type1_met": ["L1FastJet"],
            # names of the uncertainties to be applied
            "uncertainty_sources": [
                "Total",
                "CorrelationGroupMPFInSitu",
                "CorrelationGroupIntercalibration",
                "CorrelationGroupbJES",
                "CorrelationGroupFlavor",
                "CorrelationGroupUncorrelated",
            ],
        })

    :param reqs: Requirement dictionary for this
        :py:class:`~columnflow.calibration.Calibrator` instance
    :param inputs: Additional inputs, currently not used
    :param reader_targets: TODO: add documentation
    """

    bundle = reqs["external_files"]
    
    # import the correction sets from the external file
    import correctionlib

    correction_set = correctionlib.CorrectionSet.from_string(
        self.get_jec_file(bundle.files).load(formatter="gzip").decode("utf-8"),
    )

    # compute JEC keys from config information
    jec_cfg = self.get_jec_config()

    def make_jme_keys(names, jec=jec_cfg, is_data=self.dataset_inst.is_data):
        if is_data:

            jec_era = self.dataset_inst.get_aux("jec_era", None)
            # if no special JEC era is specified, infer based on 'era'
            if jec_era is None:
                jec_era = "Run" + self.dataset_inst.get_aux("era")

        return [
            f"{jec.campaign}_{jec_era}_{jec.version}_DATA_{name}_{jec.jet_type}"
            if is_data else
            f"{jec.campaign}_{jec.version}_MC_{name}_{jec.jet_type}"
            for name in names
        ]

    # take sources from constructor or config
    sources = self.uncertainty_sources
    if sources is None:
        sources = jec_cfg.uncertainty_sources
    
    if self.dataset_inst.is_data :
        jec_keys = make_jme_keys(jec_cfg.levels_DATA)
    else :
        jec_keys = make_jme_keys(jec_cfg.levels_MC)
    jec_keys_subset_type1_met = make_jme_keys(jec_cfg.levels_for_type1_met)
    junc_keys = make_jme_keys(sources, is_data=False)  # uncertainties only stored as MC keys

    # store the evaluators
    self.evaluators = {
        "jec": get_evaluators(correction_set, jec_keys),
        "jec_subset_type1_met": get_evaluators(correction_set, jec_keys_subset_type1_met),
        "junc": dict(zip(sources, get_evaluators(correction_set, junc_keys))),
    }

# custom jec calibrator that only runs nominal correction
jec_nominal = jec.derive("jec_nominal", cls_dict={"uncertainty_sources": []})
