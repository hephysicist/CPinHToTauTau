import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from httcp.util import enforce_hcand_type
ak = maybe_import("awkward")
np = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

@producer(
    uses = 
    {
        f"hcand.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {attach_coffea_behavior},
    produces={
        "hcand_obj.mass"
    },
)
def hcand_mass_(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing dilepton mass...")
    from coffea.nanoevents.methods import vector
    events = self[attach_coffea_behavior](events, **kwargs)
    lep = []
    for i in range(2):
        
        hcand_lep = events.hcand[:,i]
        lep.append( ak.zip(
            {
                "pt": hcand_lep.pt,
                "eta": hcand_lep.eta,
                "phi": hcand_lep.phi,
                "mass": hcand_lep.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        ))
    hcand_obj = lep[0] + lep[1]
    events = set_ak_column_f32(events,f"hcand_obj.mass", ak.where(hcand_obj.mass2 >=0, hcand_obj.mass, EMPTY_FLOAT))
    return events


@producer(
    uses = 
    {
        "hcand.charge",
    },
    produces={
        "hcand_obj.rel_charge"
    },
)
def rel_charge(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing  pair relative charge...")
    rel_ch = ak.prod(events.hcand.charge, axis = 1)
    events = set_ak_column_f32(events, "hcand_obj.rel_charge", rel_ch) 
    return events

@producer(
    uses = 
    {
        f"hcand.{var}" for var in ["pt","phi"]
    } | {
        f"PuppiMET.{var}" for var in ["pt","phi"] 
    } | {attach_coffea_behavior},
    produces={
        "mT"
    },
)
def mT(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("producing mT...")
    lep1 = ak.zip({
            "pt": events.hcand[:,0:1].pt,
            "eta": events.hcand[:,0:1].eta,
            "phi": events.hcand[:,0:1].phi,
            "mass": events.hcand[:,0:1].mass,
           
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=coffea.nanoevents.methods.vector.behavior,)

    cos_dphi = np.cos(lep1.delta_phi(events.PuppiMET))
    mT_values = np.sqrt(2*lep1.pt*events.PuppiMET.pt * (1 - cos_dphi))
    mT_values = ak.fill_none(mT_values, EMPTY_FLOAT)
    events = set_ak_column_f32(events,"mT", mT_values)
    return events


def hcand_mt(lep: ak.Array, MET: ak.Array) -> ak.Array:
    print("producing mT...")
    cos_dphi = np.cos(lep.delta_phi(MET))
    mT_values = np.sqrt(2*lep.pt*MET.pt * (1 - cos_dphi))
    return ak.fill_none(mT_values, EMPTY_FLOAT)

@producer(
    uses = 
    {
        f"hcand.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {
        f"Muon.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {
        f"Electron.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {
        f"Tau.{var}" for var in ["pt", "eta","phi", "mass","charge"]
    } | {attach_coffea_behavior}, # | {"channel_id"},
    produces={
        "hcand_obj.mass"
    },
)
def hcand_mass(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    print("Producing dilepton mass...")
    from coffea.nanoevents.methods import vector
    events = self[attach_coffea_behavior](events, **kwargs)
    lep = []
    for i in range(2):
        hcand_lep = events.hcand[:,i]
        lep.append( ak.zip(
            {
                "pt": hcand_lep.pt,
                "eta": hcand_lep.eta,
                "phi": hcand_lep.phi,
                "mass": hcand_lep.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
        ))
    hcand_obj = lep[0] + lep[1]
    events = set_ak_column_f32(events,f"hcand_obj.mass", ak.where(hcand_obj.mass2 >=0, hcand_obj.mass, EMPTY_FLOAT))
    return events 

    
   
