import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from MSSM_H_tt.util import get_lep_p4
ak = maybe_import("awkward")
np = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
def delta_phi(phi1, phi2):
        """Compute Δφ with proper wrapping into [−π,π]."""
        dphi = phi1 - phi2
        dphi = ak.where(dphi > np.pi,  dphi - 2*np.pi, dphi)
        dphi = ak.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
        return dphi
def mT(p4_1, p4_2):
        """Transverse mass between two four‑vectors p4_1 and p4_2."""
        dphi = delta_phi(p4_1.phi, p4_2.phi)
        return np.sqrt(2 * p4_1.pt * p4_2.pt * (1 - np.cos(dphi)))
# def hcand_mt(lep: ak.Array, MET: ak.Array) -> ak.Array:
#     print("Producing mT...")
#     delta_phi = lep.phi - MET.phi
#     delta_phi = ak.where(delta_phi > np.pi, delta_phi - 2*np.pi, delta_phi)
#     delta_phi = ak.where(delta_phi < -np.pi, delta_phi + 2*np.pi, delta_phi)
#     cos_dphi = np.cos(delta_phi)
#     mT_values = np.sqrt(2*lep.pt*MET.pt * (1 - cos_dphi))
#     return ak.fill_none(mT_values, EMPTY_FLOAT)

@producer(
    uses={
        'hcand_*', 'PuppiMET*'
    },
    produces={
        'hcand_*'
    },
)
def hcand_fields(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:
    print("producing pTll...")
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        p4 = {}
        for the_lep in hcand.fields: p4[the_lep] = get_lep_p4(hcand[the_lep]) 
        
        pair = p4['lep0'] + p4['lep1']
        dilep_mass = pair.mass
        hcand['mass'] = ak.where(dilep_mass > 0, dilep_mass , EMPTY_FLOAT)
        dilep_pt = pair.pt
        hcand['pt'] = ak.where(dilep_pt > 0, dilep_pt , EMPTY_FLOAT)
        delta_r = ak.flatten(p4['lep0'].metric_table(p4['lep1']), axis=2)
        hcand['delta_r'] = ak.where(delta_r > 0, delta_r , EMPTY_FLOAT)
        hcand['rel_charge'] = hcand.lep0.charge * hcand.lep1.charge
        
        events = set_ak_column(events, f'hcand_{ch_str}', hcand) 
    return events

@producer(
    uses={
        'hcand_*', 'PuppiMET*'
    },
    produces={
        'hcand_*'
    },
)
def hcand_mt(self: Producer, 
             events: ak.Array,
             **kwargs
             ) -> ak.Array:
    print("producing mT...")
    channels = self.config_inst.channels.names()
    ch_objects = self.config_inst.x.ch_objects
    MET = events.PuppiMET  
    for ch_str in channels:
        hcand = events[f'hcand_{ch_str}']
        p4 = {}
        # met_filled = ak.fill_none(events.PuppiMET, 0)
        for the_lep in hcand.fields:
            if the_lep in ['lep0', 'lep1']:
                p4[the_lep] = get_lep_p4(hcand[the_lep])
        if ch_str in ['etau', 'mutau']:
            mt = mT(p4['lep0'], MET)
            hcand['mt'] = ak.where(mt > 0, mt , EMPTY_FLOAT)
        if ch_str == 'emu':
            mt_e   = mT(p4['lep0'],MET)  
            mt_mu  = mT(p4['lep1'],MET)           
            mt_emu = mT(p4['lep0'],p4['lep1'])  
            mt_tot = np.sqrt(mt_e**2 + mt_mu**2 + mt_emu**2)
            hcand['mt_e'] = ak.where(mt_e > 0, mt_e , EMPTY_FLOAT)
            hcand['mt_mu'] = ak.where(mt_mu > 0, mt_mu , EMPTY_FLOAT)
            hcand['mt_emu'] = ak.where(mt_emu > 0, mt_emu , EMPTY_FLOAT)
            hcand['mt_tot'] = ak.where(mt_tot > 0, mt_tot , EMPTY_FLOAT)
        
        events = set_ak_column(events, f'hcand_{ch_str}', hcand) 

    return events 

    
   
