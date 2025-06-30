import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column, EMPTY_FLOAT, Route, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from httcp.util import get_lep_p4
ak = maybe_import("awkward")
np = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)

def calc_mt(obj1: ak.Array, obj2: ak.Array)-> ak.Array:
        d_phi = obj1.phi - obj2.phi
        d_phi = ak.where(abs(d_phi) > np.pi,
                             d_phi -2 * np.sign(d_phi) * np.pi,
                             d_phi)
        cos_dphi = np.cos(d_phi)
        mt_values = np.sqrt(2*obj1.pt*obj2.pt * (1 - cos_dphi))
        mt_values = ak.nan_to_num(mt_values,EMPTY_FLOAT)
        return mt_values
    
def calc_pt(obj1: ak.Array, obj2: ak.Array)-> ak.Array:
        d_phi = obj1.phi - obj2.phi
        d_phi = ak.where(abs(d_phi) > np.pi,
                             d_phi -2 * np.sign(d_phi) * np.pi,
                             d_phi)
        cos_dphi = np.cos(d_phi)
        pt_values = np.sqrt(obj1.pt**2 + obj2.pt**2 + 2 * obj1.pt * obj2.pt * cos_dphi)
        pt_values = ak.nan_to_num(pt_values, EMPTY_FLOAT)
        return pt_values

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
    ch_str = self.config_inst.channels.names()[0]
    hcand = events[f'hcand_{ch_str}']
    lep0_p4 = get_lep_p4(hcand.lep0)
    lep1_p4 = get_lep_p4(hcand.lep1)
    pair_p4 = lep0_p4+lep1_p4
    met = events.PuppiMET
    
    hcand['mass']   = ak.where(pair_p4.mass > 0, pair_p4.mass , EMPTY_FLOAT)
    hcand['pt_vis'] = ak.where(pair_p4.pt > 0, pair_p4.pt , EMPTY_FLOAT)
    hcand['pt_tt']  = calc_pt(pair_p4, met)
    
    delta_r = ak.flatten(lep0_p4.metric_table(lep1_p4), axis=2)
    hcand['delta_r'] = ak.where(delta_r > 0, delta_r , EMPTY_FLOAT)
    
    delta_eta = lep0_p4.eta - lep1_p4.eta
    hcand['delta_eta']  = ak.where(abs(delta_eta) < 10, delta_eta , EMPTY_FLOAT)
    hcand['delta_phi']  = lep0_p4.deltaphi(lep1_p4)
    
    hcand['rel_charge'] = hcand.lep0.charge * hcand.lep1.charge
    
    hcand['mt_0']    = calc_mt(lep0_p4, met)
    hcand['mt_1']    = calc_mt(lep1_p4, met)
    hcand['mt_ll']   = calc_mt(pair_p4, met)
    hcand['mt_vis']  = calc_mt(lep0_p4,lep1_p4) 
    
    
    mt_tot_mask = (hcand['mt_0']>0) & (hcand['mt_1']>0) & (hcand['mt_vis']>0)
    hcand['mt_tot'] = ak.where(mt_tot_mask,
                        np.sqrt(hcand['mt_0']**2 + hcand['mt_1']**2 + hcand['mt_vis']**2),
                        EMPTY_FLOAT)
    hcand['delta_phi_0_met'] = lep0_p4.deltaphi(met)
    hcand['delta_phi_1_met'] = lep1_p4.deltaphi(met)
    events = set_ak_column(events, f'hcand_{ch_str}', hcand) 
    return events