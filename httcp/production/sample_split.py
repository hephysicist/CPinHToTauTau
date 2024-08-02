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
    uses={f"Tau.{var}" for var in [
                "decayMode", "genPartFlav"
                ] 
    } | {"process_id"},
    produces={
        "process_id"
    },
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
  
    dm = flat_np_view(events.Tau.decayMode, axis=1)
    match = flat_np_view(events.Tau.genPartFlav, axis=1)
    process_id = np.array(events.process_id)
    
    channles_mask = {}
    for the_ch in ["etau","mutau", "tautau"]: channles_mask[the_ch] = events.channel_id == self.config_inst.get_channel(the_ch).id
   
    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    
    mu2tau_fakes_mask = (match == tau_part_flav["prompt_mu"])
    mu_tau_fake_id = 51001*ak.ones_like(events.process_id) 
    mu_tau_gen_id = 51002*ak.ones_like(events.process_id) 
    process_id[mu2tau_fakes_mask] = mu_tau_fake_id[mu2tau_fakes_mask]
    process_id[~mu2tau_fakes_mask] = mu_tau_gen_id[~mu2tau_fakes_mask]
    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int32)
    return events


