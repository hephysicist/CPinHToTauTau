# coding: utf-8

"""
Trigger selection methods.
"""

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, optional_column as opt

np = maybe_import("numpy")
ak = maybe_import("awkward")


@selector(
    uses={
        # nano columns
        "TrigObj.id", "TrigObj.pt", "TrigObj.eta", "TrigObj.phi", "TrigObj.filterBits","HLT.*"
    },
    produces={
        # new columns
        "trigger_ids",
    },
    exposed=True,
)
def trigger_selection(
    self: Selector,
    events: ak.Array,
    **kwargs,
) -> tuple[ak.Array, SelectionResult]:
    """
    HLT trigger path selection.
    """
    any_fired = False
    any_fired_all_legs_match = False
    trigger_data = []
    trigger_ids = []

    # index of TrigObj's to repeatedly convert masks to indices
    index = ak.local_index(events.TrigObj)

    for trigger in self.config_inst.x.triggers:
        # skip the trigger if it does not apply to the dataset
        #if not trigger.applies_to_dataset(self.dataset_inst):
        #    continue
        #perform leg-specific selections first
        leg_masks = []
        all_legs_match = True
        
        for leg in trigger.legs:
            # start with a True mask
            leg_mask = abs(events.TrigObj.id) >= 0
            leg_mask = ak.enforce_type(leg_mask, f"var * bool")
            # pdg id selection
            if leg.pdg_id is not None:
                Trig_pdg_ID = (abs(events.TrigObj.id) == leg.pdg_id)
                Trig_pdg_ID = ak.enforce_type(Trig_pdg_ID, f"var * bool")
                leg_mask = leg_mask & Trig_pdg_ID
            # # pt cut
            # if leg.min_pt is not None:
            #     Trig_pt = (events.TrigObj.pt >= leg.min_pt)
            #     Trig_pt = ak.enforce_type(Trig_pt, f"var * bool")
            #     leg_mask = leg_mask & (events.TrigObj.pt >= leg.min_pt)
            # # eta cut
            # if leg.min_eta is not None:
            #     Trig_eta = (abs(events.TrigObj.eta) < leg.min_eta)
            #     Trig_eta = ak.enforce_type(Trig_eta, f"var * bool")
            #     leg_mask = leg_mask & Trig_eta
            # trigger bits match
            if leg.trigger_bits is not None:
                bit_mask = 0
                for bit in leg.trigger_bits:
                    bit_mask += 1<<(bit-1)
                    Filter_Bits = ((events.TrigObj.filterBits & bit_mask) != 0)
                    Filter_Bits =  ak.enforce_type(Filter_Bits, f"var * bool")
                    leg_mask = leg_mask & Filter_Bits
            leg_masks.append(index[leg_mask])
            # at least one object must match this leg
            all_legs_match = all_legs_match & ak.any(leg_mask, axis=1)
        
       
        # get bare decisions
        fired = events.HLT[trigger.hlt_field]
        any_fired = any_fired | fired
        # final trigger decision
        fired_and_all_legs_match = fired & all_legs_match
        any_fired_all_legs_match = (any_fired_all_legs_match | fired_and_all_legs_match)
        # store all intermediate results for subsequent selectors
        trigger_data.append((trigger, fired_and_all_legs_match, leg_masks))

        # store the trigger id
        ids = ak.where(fired_and_all_legs_match, np.float64(trigger.id), np.float64(np.nan))
        trigger_ids.append(ak.singletons(ak.nan_to_none(ids)))

    # store the fired trigger ids
    trigger_ids = ak.concatenate(trigger_ids, axis=1)
    events = set_ak_column(events, "trigger_ids", trigger_ids, value_type=np.int64)
    return events, SelectionResult(
        steps={
            "trigger": any_fired,
        },
        aux={
            "trigger_data": trigger_data,
        },
    )


@trigger_selection.init
def trigger_selection_init(self: Selector) -> None:
    if getattr(self, "dataset_inst", None) is None:
        return

    # full used columns
    self.uses |= {
        opt(trigger.name)
        for trigger in self.config_inst.x.triggers
        if trigger.applies_to_dataset(self.dataset_inst)
    }




