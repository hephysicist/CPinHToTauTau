"""
Column production methods related to higher-level features.
"""
import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional
from MSSM_H_tt.util import get_lep_p4, get_ip_p4
from columnflow.production.util import attach_coffea_behavior

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")
warn = maybe_import("warnings")

###################
## Generator Z pt
###################
@producer(
    uses={
        "GenPart.*"
    },
    produces={
        "GenZ.pt", "GenZ.eta", "GenZ.phi", "GenZ.mass"
    },
    mc_only=True,
)
def generatorZ(
        self : Producer,
        events : ak.Array,
        **kwargs) -> ak.Array:


    genpart_indices = ak.local_index(events.GenPart.pt)

    sel_gen_ids = genpart_indices[
        ((((np.abs(events.GenPart.pdgId) >= 11) & (np.abs(events.GenPart.pdgId) <= 16)) & (events.GenPart.statusFlags & 256 != 0)) & (events.GenPart.status == 1)) | (events.GenPart.statusFlags & 1024 != 0)]

    gen_part = events.GenPart[sel_gen_ids]

    # form LV
    gen_part = ak.Array(gen_part, behavior=coffea.nanoevents.methods.nanoaod.behavior)
    p4_gen_part = ak.with_name(gen_part, "PtEtaPhiMLorentzVector")    
    p4_gen_part = ak.zip(
        {
            "x" : p4_gen_part.px,
            "y" : p4_gen_part.py,
            "z" : p4_gen_part.pz,
            "t" : p4_gen_part.energy,
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )

    sum_p4_gen_part = ak.zip(
        {
            "x" : ak.sum(p4_gen_part.x, axis=1),
            "y" : ak.sum(p4_gen_part.y, axis=1),
            "z" : ak.sum(p4_gen_part.z, axis=1),
            "t" : ak.sum(p4_gen_part.energy, axis=1),
        },
        with_name="LorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,
    )

    p4_gen_part_array = ak.zip(
        {
            "pt"  : ak.nan_to_num(sum_p4_gen_part.pt, 0.0),  # nan values are for empty genpart
            "eta" : ak.nan_to_num(sum_p4_gen_part.eta, 0.0), 
            "phi" : ak.nan_to_num(sum_p4_gen_part.phi, 0.0),
            "mass": ak.nan_to_num(sum_p4_gen_part.mass, 0.0),
        }
    )

    events = set_ak_column(events, "GenZ", p4_gen_part_array)
    return events
