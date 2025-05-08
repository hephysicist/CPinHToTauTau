import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={
        "Jet*", "Jet.hadronFlavour",
        },
    produces={
         f"btag_weight_SF_{shift}"
        for shift in ["nom"]
    },
    mc_only=True,
)
def btag_weight_SF(
    self: Producer, 
    events: ak.Array, 
    do_syst: bool,
    **kwargs,  
) -> ak.Array:
    
    shifts = ["central"]
    if do_syst: shifts=[*shifts] 
    sf_values = {}
    tag = self.config_inst.x.tag
    year = self.config_inst.x.year

    jet_mask = ((events.Jet.pt >= 20) & 
                (abs(events.Jet.eta) < 2.5) & 
                (events.Jet.jetId & 0b10 == 0b10))
    
    #Removing NaNs from discriminat
    dis = events.Jet.btagDeepFlavB 
    nan_mask = np.isnan(dis)
    mask = ~np.isnan(dis)
    jet_mask = jet_mask & mask
    
    Jet = events.Jet[jet_mask]

    for the_shift in shifts: sf_values[the_shift] = np.ones_like(events.event, dtype=np.float32)
    # Create sf array template to make copies and dict for finnal results of all systematics

    flavor = Jet.hadronFlavour
    eta = abs(Jet.eta)
    pt = Jet.pt
    discriminant =Jet.btagDeepFlavB

    #Prepare a tuple with the inputs of the correction evaluator
    btag_sf_args = lambda syst : (syst,flavor,eta,pt,discriminant)

    #Loop over the shifts and calculate for each shift btag scale factor
    for the_shift in shifts:
        sf = ak.ones_like(pt)
        sf = sf * self.btag_sf_corr.evaluate(*btag_sf_args(the_shift))
        sf_values[the_shift] = sf_values[the_shift] * ak.fill_none(sf, 1.)
                    
    rename_systs = {"central" : "nom",} 

    selection_stats = self.task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    ) 
    pid = self.dataset_inst.processes.ids()[0]
    pid = str(pid)
    sum_mc_weight_before_btag_weigth = selection_stats["sum_mc_weight_selected_per_process"][pid]
    sum_mc_weight_after_btag_weigth = selection_stats["sum_mc_weight_selected_with_btag_weight_per_process"][pid]
    r = sum_mc_weight_before_btag_weigth/sum_mc_weight_after_btag_weigth

    for the_shift in shifts: 
        btag_w = sf_values[the_shift]
        w_event = ak.prod(btag_w,axis=1)

        events = set_ak_column_f32(events, f"btag_weight_SF_{rename_systs[the_shift]}", w_event*r)
        
    return events

@btag_weight_SF.requires
def btag_weight_SF_requires(self: Producer, reqs: dict) -> None:
    from columnflow.tasks.selection import MergeSelectionStats
    reqs["selection_stats"] = MergeSelectionStats.req_different_branching(
        self.task,
        branch=-1 if self.task.is_workflow() else 0,
    )
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@btag_weight_SF.setup
def btag_weight_SF_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    
    selection_stats = self.task.cached_value(
        key="selection_stats",
        func=lambda: inputs["selection_stats"]["stats"].load(formatter="json"),
    )

    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.btag_sf_corr.load(formatter="gzip").decode("utf-8"),
    )
    self.btag_sf_corr = correction_set[self.config_inst.x.btag_sf[0]]

    # nan_mask = np.isnan(discriminant)
    
    # mask = ~np.isnan(discriminant)

    # # nan_mask = np.isnan(discriminant)
    # if ak.any(nan_mask):
    #     print("NaN values found in the array.")
    # else:
    #     print("No NaN values in the array.")
    # nan_count = ak.sum(np.isnan(discriminant))
    # print(f"Found {nan_count} NaN values.")
    
    # # Filter out the nested subarrays that contain NaN values.
    # filtered = discriminant[mask]

    # print("Original array:")
    # print(discriminant)
    # print("\nMask for subarrays without NaNs:")
    # print(mask)
    # print("\nFiltered array (NaN-containing entries removed):")
    # print(filtered)

    # flavor_ = flavor[mask]
    # eta_ = eta[mask]
    # pt_ = pt[mask]
    # discriminant_ = discriminant[mask]