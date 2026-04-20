import functools
from columnflow.production import Producer, producer
from columnflow.columnar_util import set_ak_column, has_ak_column, flat_np_view, optional_column as optional
from columnflow.util import maybe_import, load_correction_set

import law
ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses={'event', optional("LHEScaleWeight.*"),optional("PSWeight.*")},
    produces={"lhe_weight*", "ps_weight*"},
    mc_only=True,
)
def theor_unc(self: Producer, events: ak.Array, **kwargs) -> ak.Array:

    
    lhe_weight = {}
    lhe_weight['nominal'] = ak.ones_like(events.event, dtype=np.float32)        
    unit_weight = ak.ones_like(events.event, dtype=np.float32)   
    events = set_ak_column(events, "lhe_weight", unit_weight)
    for syst_name, bin in self.config_inst.x.lhe_variations.items():
        if "LHEScaleWeight" in events.fields: 
            mask = ak.num(events.LHEScaleWeight, axis = 1) == 9 #Total number of LHE weights per events, if less, then don't apply any weights 
            the_weight = ak.fill_none(ak.firsts(events.LHEScaleWeight[:,bin:]),1)
            lhe_weight[syst_name] = ak.where(mask,
                                                the_weight,
                                                unit_weight)
        else:
            lhe_weight[syst_name] = unit_weight
        events = set_ak_column(events, f"lhe_weight_{syst_name}", lhe_weight[syst_name])
                
    ps_weight = {}
    ps_weight['nominal'] = ak.ones_like(events.event, dtype=np.float32) 
    events = set_ak_column(events, "ps_weight", unit_weight)
    for syst_name, bin in self.config_inst.x.ps_variations.items():
        if "PSWeight" in events.fields:
            mask = ak.num(events.PSWeight, axis = 1) == 4 #Total number of PS weights per events, if less, then don't apply any weights 
            the_weight = ak.fill_none(ak.firsts(events.PSWeight[:,bin:]),1)
            ps_weight[syst_name] = ak.where(mask,
                                                the_weight,
                                                unit_weight)
        else:
            ps_weight[syst_name] = unit_weight
        events = set_ak_column(events, f"ps_weight_{syst_name}", ps_weight[syst_name])
    return events










# def Signal_Theory(events, weights, year, dataset_name, is_correction=False, **kwargs):

#     if not is_correction:

#         all_variations = [
#             "Scale_muR",
#             "Scale_muF",
#             "PS_ISR",
#             "PS_FSR"
#         ]

#         variations = {}
#         for syst_name in all_variations:
#             variations[syst_name] = {}
#             variations[syst_name]["up"] = np.ones(len(weights._weight))
#             variations[syst_name]["down"] = np.ones(len(weights._weight))

#         if "LHEScaleWeight" in events.fields:
#             LHEScaleWeight = events.LHEScaleWeight

#             # pad to 9 and fill missing with 1.0; clip avoids extending lists longer than 9
#             padded = ak.pad_none(LHEScaleWeight, 9, axis=1, clip=True)

#             counts = ak.num(LHEScaleWeight, axis=1)
#             bad_mask = counts != 9
#             if ak.any(bad_mask):
#                 logger.warning(f"{ak.sum(bad_mask)} events have != 9 LHEScaleWeight entries; missing indices set to 1.0")

#             scale_muR_0p5 = ak.fill_none(padded[:, 1], 1.0)  # muR 0.5, muF 1.0
#             scale_muR_2p0 = ak.fill_none(padded[:, 7], 1.0)  # muR 2.0, muF 1.0
#             scale_muF_0p5 = ak.fill_none(padded[:, 3], 1.0)  # muR 1.0, muF 0.5
#             scale_muF_2p0 = ak.fill_none(padded[:, 5], 1.0)  # muR 1.0, muF 2.0
#         else:
#             logger.warning("LHEScaleWeight not found; setting all scale variations to 1.0")
#             scale_muR_0p5 = np.ones(len(weights._weight))
#             scale_muR_2p0 = np.ones(len(weights._weight))
#             scale_muF_0p5 = np.ones(len(weights._weight))
#             scale_muF_2p0 = np.ones(len(weights._weight))

#         variations["Scale_muR"]["up"] *= np.asarray(scale_muR_2p0)
#         variations["Scale_muR"]["down"] *= np.asarray(scale_muR_0p5)
#         variations["Scale_muF"]["up"] *= np.asarray(scale_muF_2p0)
#         variations["Scale_muF"]["down"] *= np.asarray(scale_muF_0p5)

#         if "PSWeight" in events.fields:
#             PSWeight = events.PSWeight

#             # pad to 4 slots and fill with 1.0 where missing
#             ps_padded = ak.pad_none(PSWeight, 4, axis=1, clip=True)

#             ps_counts = ak.num(PSWeight, axis=1)
#             ps_bad_mask = ps_counts < 4
#             if ak.any(ps_bad_mask):
#                 logger.warning(f"{ak.sum(ps_bad_mask)} events have < 4 PSWeight entries; missing indices set to 1.0")

#             ps_ISR_2p0 = ak.fill_none(ps_padded[:, 0], 1.0)  # ISR = 2, FSR = 1
#             ps_FSR_2p0 = ak.fill_none(ps_padded[:, 1], 1.0)  # ISR = 1, FSR = 2
#             ps_ISR_0p5 = ak.fill_none(ps_padded[:, 2], 1.0)  # ISR = 0.5, FSR = 1
#             ps_FSR_0p5 = ak.fill_none(ps_padded[:, 3], 1.0)  # ISR = 1, FSR = 0.5
#         else:
#             logger.warning("PSWeight not found; setting all parton-shower variations to 1.0")
#             ps_ISR_2p0 = np.ones(len(weights._weight))
#             ps_ISR_0p5 = np.ones(len(weights._weight))
#             ps_FSR_2p0 = np.ones(len(weights._weight))
#             ps_FSR_0p5 = np.ones(len(weights._weight))

#         variations["PS_ISR"]["up"] *= np.asarray(ps_ISR_2p0)
#         variations["PS_ISR"]["down"] *= np.asarray(ps_ISR_0p5)
#         variations["PS_FSR"]["up"] *= np.asarray(ps_FSR_2p0)
#         variations["PS_FSR"]["down"] *= np.asarray(ps_FSR_0p5)

#         sfs_up = [variations[syst_name]["up"] for syst_name in all_variations]
#         sfs_down = [variations[syst_name]["down"] for syst_name in all_variations]

#         weights.add_multivariation(
#             name="Signal_Theory",
#             weight=np.ones(len(weights._weight)),
#             weightsUp=sfs_up,
#             weightsDown=sfs_down,
#             modifierNames=all_variations
#         )

#         return weights
