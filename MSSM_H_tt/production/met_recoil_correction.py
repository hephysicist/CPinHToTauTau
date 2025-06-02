import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask

from httcp.util import get_lep_p4

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)



@producer(
    uses={
        "GenPart.*"
    },
    produces=
    {
        f"{the_part}*" 
        for the_part in ['gen_boson','gen_boson_vis']
    },
    mc_only=True,
)

def gen_boson(
        self : Producer,
        events : ak.Array,
        **kwargs) -> ak.Array:

    gen_part = events.GenPart
    genpart_indices = ak.local_index(gen_part.pt)

    pdg_id = abs(gen_part.pdgId)
    
    #helper function to check a specific bit in the number
    check_bit = lambda num, bit_pos: num>>bit_pos & 1
    # get all decay products of leptonic W, Z, or H decays
    
    #Choose electrons (11) or muons (13) that marked as stable (status == 1) comming fromHardProcess (bit 8)
    is_e_or_mu = (((pdg_id == 11) | (pdg_id == 13)) & (gen_part.status == 1)) & check_bit(gen_part.statusFlags,8)
    
    #Choose neutrinos (12,14,16) that are marked as stable (status == 1) 
    is_nu = ((pdg_id == 12) | (pdg_id == 14)| (pdg_id == 16)) & (gen_part.status == 1)
    
    #Choose tauh decay products (pions and neutrinos): isDirectHardProcessTauDecayProduct (bit 10)      
    is_tauh_decay_prod = check_bit(gen_part.statusFlags,10) 
    from_hard_proc =  check_bit(gen_part.statusFlags,8) 
    #Choose all leptons and neutrinos for the full boson momentum   
    full_mask = ((is_e_or_mu | is_nu) & from_hard_proc) | is_tauh_decay_prod
    vis_mask =  ((is_e_or_mu & from_hard_proc) | is_tauh_decay_prod) & (~is_nu)
    zeros_array =  ak.zeros_like(gen_part)
    #Select full gen objects set and visible objects set
    gen_obj = get_lep_p4(ak.where(full_mask, gen_part, zeros_array))
    vis_obj = get_lep_p4(ak.where(vis_mask, gen_part, zeros_array))
    
    events = set_ak_column(events, "gen_boson", gen_obj)
    events = set_ak_column(events, "gen_boson_vis", vis_obj)
    
    return events

### MET RECOIL CORRECTION CALCULATOR ###

@producer(
    uses={
        'PuppiMET*','hcand*', 'n_jets', 'event', optional('gen_boson*')
    },
    produces={
        'PuppiMET*'
    },
)
def met_recoil(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    met_recoil_dict = self.config_inst.x.met_recoil.datasets
    dataset_name = self.dataset_inst.name
    if dataset_name in met_recoil_dict.keys():
        print(f"Applying MET recoil correction for {dataset_name}")
        boson_p4     = events.gen_boson.sum(axis=1)
        boson_p4_vis = events.gen_boson_vis.sum(axis=1)
        zeros_arr = ak.zeros_like(events.PuppiMET.pt)
        met_p4 = ak.zip(
            {
                "pt"  : events.PuppiMET.pt,
                "eta" : zeros_arr,
                "phi" : events.PuppiMET.phi,
                "mass": zeros_arr,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=coffea.nanoevents.methods.vector.behavior,
        )
        u = met_p4 + boson_p4_vis - boson_p4
        dphi_v_u = u.phi - boson_p4.phi
        u_para = u.pt * np.cos(dphi_v_u)
        u_perp = u.pt * np.sin(dphi_v_u)
        
        ch_str = self.config_inst.channels.names()[0]
        hcand = events[f'hcand_{ch_str}']
        
        order      = met_recoil_dict[dataset_name]
        njet       = events.n_jets 
        ptll       = flat_np_view(hcand.pt)
        
        u_para_args = lambda : (order,
                                njet,
                                ptll,
                                "Upara",
                                u_para)
        
        u_perp_args = lambda : (order,
                                njet,
                                ptll,
                                "Uperp",
                                u_perp)
        
        u_para_corr = self.QuantileHistCorr.evaluate(*u_para_args())
        u_perp_corr = self.QuantileHistCorr.evaluate(*u_perp_args())
        
        dphi_v_u_corr = np.atan2(u_perp_corr, u_para_corr) 
        u_corr = ak.zip(
            {
                "pt": np.sqrt(u_para_corr**2 + u_perp_corr**2),
                "eta": zeros_arr,
                "phi": dphi_v_u_corr + boson_p4.phi,
                "mass": zeros_arr,
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        met_p4_corr = u_corr - boson_p4_vis + boson_p4  
        
        events = set_ak_column_f32(events, "PuppiMET.pt_no_recoil", events.PuppiMET.pt)
        events = set_ak_column_f32(events, "PuppiMET.phi_no_recoil", events.PuppiMET.phi)
        events = set_ak_column_f32(events, "PuppiMET.pt_recoil_effect", met_p4_corr.pt - events.PuppiMET.pt)
        events = set_ak_column_f32(events, "PuppiMET.phi_recoil_effect", met_p4_corr.phi - events.PuppiMET.phi)
        events = set_ak_column_f32(events, "PuppiMET.pt", met_p4_corr.pt)
        events = set_ak_column_f32(events, "PuppiMET.phi", met_p4_corr.phi)
    else:
        print(f"NOT applying MET recoil correction for {dataset_name}")
        events = set_ak_column_f32(events, "PuppiMET.pt_no_recoil", events.PuppiMET.pt)
        events = set_ak_column_f32(events, "PuppiMET.phi_no_recoil", events.PuppiMET.phi)
        events = set_ak_column_f32(events, "PuppiMET.pt_recoil_effect", ak.zeros_like(events.PuppiMET.pt))
        events = set_ak_column_f32(events, "PuppiMET.phi_recoil_effect", ak.zeros_like(events.PuppiMET.phi))
    
    nan_mask = ak.zeros_like(events.PuppiMET.pt, dtype=np.bool_)
    for the_name in events.PuppiMET.fields:
        the_var = events.PuppiMET[the_name]
        nan_mask = nan_mask | np.isnan(the_var)
    if ak.any(nan_mask):
        # print(f'Found nans in events.PuppiMET object:')
        # for the_evt in events[nan_mask]:
        #     print(the_evt)
        masked_met = events.PuppiMET
        for the_var in masked_met.fields:
            masked_met[the_var] = ak.where(nan_mask, EMPTY_FLOAT, masked_met[the_var])
        events = set_ak_column_f32(events, "PuppiMET", masked_met)
    return events

@met_recoil.requires
def met_recoil_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@met_recoil.setup
def met_recoil_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.met_recoil.load(formatter="gzip").decode("utf-8"),
    )
    self.QuantileHistCorr = correction_set["Recoil_correction_QuantileMapHist"]