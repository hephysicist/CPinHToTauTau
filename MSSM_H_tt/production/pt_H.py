"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict
from MSSM_H_tt.util import get_lep_p4, find_fields_with_nan, get_p2
from columnflow.columnar_util import EMPTY_FLOAT

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
import functools
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)

@producer(
    uses={ f"PuppiMET.{var}" for var in 
        ["pt", "phi",]} | {"hcand_*"},
    produces={"pt_H"},
    exposed=False,
)
def pt_H(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces 'pt_H' column to be used in categorisation
    """

    electron    = events.hcand_emu.lep0
    muon        = events.hcand_emu.lep1
    
    somma = get_p2(electron)+get_p2(muon)+get_p2(events.PuppiMET)
    pt_H = somma.rho

    events = set_ak_column(events, "pt_H", pt_H)
    return events
