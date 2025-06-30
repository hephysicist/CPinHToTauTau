import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column,remove_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses={'hcand_*', 'event'
    } | {"process_id"},
    produces={
        "process_id"
    },
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    channels =  self.config_inst.channels.names()
    process_id = events.process_id
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        tau_dm = ak.firsts(hcand.lep1.genPartFlav)
        fake_mask = ak.zeros_like(tau_dm, dtype=np.bool_)
        if ch_str == 'etau':
            fake_mask = fake_mask | (tau_dm == tau_part_flav["prompt_e"])
            fake_mask = fake_mask | (tau_dm == tau_part_flav["tau->e"])
            process_id = ak.where(ak.fill_none(fake_mask, False), 51001, process_id) #z->ee events
        elif ch_str == 'mutau':
            fake_mask = fake_mask | (tau_dm == tau_part_flav["prompt_mu"])
            fake_mask = fake_mask | (tau_dm == tau_part_flav["tau->mu"])
            process_id = ak.where(ak.fill_none(fake_mask, False), 51004, process_id) #z->mumu events
        
        process_id = ak.where(ak.fill_none(tau_dm == tau_part_flav["tau_had"], False), 51005, process_id) #z->tautau events    
    
    events = set_ak_column(events, "process_id", process_id, value_type=np.int64)
    return events


