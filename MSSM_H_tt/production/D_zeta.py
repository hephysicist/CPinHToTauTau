"""
Produce channel_id column. This function is called in the main selector
"""

from columnflow.production import Producer, producer
from columnflow.selection import Selector, SelectionResult, selector
from columnflow.columnar_util import set_ak_column
from columnflow.util import maybe_import, DotDict
from MSSM_H_tt.util import get_lep_p4, find_fields_with_nan, get_p2

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
    produces={"D_zeta"},
    exposed=False,
)
def D_zeta(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> ak.Array:
    """
    This function that produces 'D_zeta' column to be used in categorisation
    """

    electron    = events.hcand_emu.lep0
    muon        = events.hcand_emu.lep1
    
    puppi_met_p2 = get_p2(events.PuppiMET)  
    electron_p2 = get_p2(electron)
    muon_p2 = get_p2(muon)
    
    sum_electron_muon = electron_p2 + muon_p2
    
    zeta =  sum_electron_muon/sum_electron_muon.rho

    p_zeta_vis = sum_electron_muon.dot(zeta)
    
    p_zeta_miss = puppi_met_p2.dot(zeta)
    
    d_zeta = p_zeta_miss - 0.85*p_zeta_vis

    events = set_ak_column(events, "D_zeta", d_zeta)
    return events
