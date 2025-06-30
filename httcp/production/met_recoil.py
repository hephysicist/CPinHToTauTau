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


#    "name": "Recoil_correction_QuantileMapHist",
#    "description": "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_QuantileMapHist' returns Upara/Uperp using the CDF obtained directly from histogram instead.",
#    "version": 1,
#    "inputs": [
#     {
#      "name": "order",
#      "type": "string",
#      "description": "Order of samples: LO, NLO, NNLO"
#     },
#     {
#      "name": "njet",
#      "type": "real",
#      "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"
#     },
#     {
#      "name": "ptll",
#      "type": "real",
#      "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"
#     },
#     {
#      "name": "var",
#      "type": "string",
#      "description": "The variable name you are giving the input for: 'Upara', 'Uperp' (string). The output will be for the same kind of variable."
#     },
#     {
#      "name": "val",
#      "type": "real",
#      "description": "Input value of either Upara or Uperp"
#     }



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
    
    

# def correct_met_recoil(met: awkward.Array, pairs: awkward.Array, era: str, isLO: bool) -> awkward.Array:

#     met_pt = met.pt
#     met_phi = met.phi
#     gen_boson_pt = pairs.gen_boson_pT
#     gen_boson_phi = pairs.gen_boson_phi
#     gen_boson_visible_pt = pairs.gen_boson_visible_pT
#     gen_boson_visible_phi = pairs.gen_boson_visible_phi

#     met_x = met_pt * np.cos(met_phi)
#     met_y = met_pt * np.sin(met_phi)
#     gen_boson_x = gen_boson_pt * np.cos(gen_boson_phi)
#     gen_boson_y = gen_boson_pt * np.sin(gen_boson_phi)
#     gen_boson_visible_x = gen_boson_visible_pt * np.cos(gen_boson_visible_phi)
#     gen_boson_visible_y = gen_boson_visible_pt * np.sin(gen_boson_visible_phi)

#     U_x = met_x + gen_boson_visible_x - gen_boson_x
#     U_y = met_y + gen_boson_visible_y - gen_boson_y

#     U_pt = np.sqrt(U_x**2 + U_y**2)
#     U_phi = np.arctan2(U_y, U_x)

#     Upara = U_pt * np.cos(U_phi - gen_boson_phi)
#     Uperp = U_pt * np.sin(U_phi - gen_boson_phi)

#     higgs_dna_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
#     path_to_json = os.path.join(higgs_dna_path, "systematics/ditau/JSONs/MET_Recoil/Recoil_corrections.json.gz")
#     evaluator = correctionlib.CorrectionSet.from_file(path_to_json)["Recoil_correction_QuantileMapHist"]

#     n_jets = awkward.to_numpy(pairs.n_jets_recoil).astype(float)

#     if isLO:
#         Upara_corr = evaluator.evaluate(era, "LO", n_jets, pairs.gen_boson_pT, "Upara", Upara)
#         Uperp_corr = evaluator.evaluate(era, "LO", n_jets, pairs.gen_boson_pT, "Uperp", Uperp)
#     else:
#         Upara_corr = evaluator.evaluate(era, "NLO", n_jets, pairs.gen_boson_pT, "Upara", Upara)
#         Uperp_corr = evaluator.evaluate(era, "NLO", n_jets, pairs.gen_boson_pT, "Uperp", Uperp)

#     U_pt_corr = np.sqrt(Upara_corr**2 + Uperp_corr**2)
#     U_phi_corr = np.arctan2(Uperp_corr, Upara_corr) + gen_boson_phi
#     U_x_corr = U_pt_corr * np.cos(U_phi_corr)
#     U_y_corr = U_pt_corr * np.sin(U_phi_corr)
#     met_x_corr = U_x_corr - gen_boson_visible_x + gen_boson_x
#     met_y_corr = U_y_corr - gen_boson_visible_y + gen_boson_y

#     met_pt_corr = np.sqrt(met_x_corr**2 + met_y_corr**2)
#     met_phi_corr = np.arctan2(met_y_corr, met_x_corr)

#     met["pt"] = met_pt_corr
#     met["phi"] = met_phi_corr

#     return met























# # "name": "Recoil_correction_Uncertainty",
# #    "description": "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_Uncertainty' gives uncertainties which are to be applies in addition to the nominal corrections. Note that a different input is required. The uncertainties were derived based on corrections applied with the QuantileMapHist method.",
# #    "version": 1,
# #    "inputs": [
# #     {
# #      "name": "order",
# #      "type": "string",
# #      "description": "Order of samples: LO, NLO, NNLO"
# #     },
# #     {
# #      "name": "njet",
# #      "type": "real",
# #      "description": "Number of jets with pT>30 GeV at |eta|<2.5, plus jets with pT>50 GeV outside of tracker region (must be converted to float value for technical reasons)"
# #     },
# #     {
# #      "name": "ptll",
# #      "type": "real",
# #      "description": "Full gen-level boson pT, obtained from all gen-level decay products, including neutrinos"
# #     },
# #     {
# #      "name": "var",
# #      "type": "string",
# #      "description": "The variable name you are giving the input for: 'Hpara', 'Hperp' (string). Note that this is a DIFFERENT vairable from Upara/Uperp! The output will be for the same kind of variable."
# #     },
# #     {
# #      "name": "val",
# #      "type": "real",
# #      "description": "Input value of either Hpara or Hperp"
# #     },
# #     {
# #      "name": "syst",
# #      "type": "string",
# #      "description": "Systematic variations on Response and Resolution: 'RespUp', 'RespDown', 'ResolUp', 'ResolDown'"
# #     }


# ### DY PTLL RECOIL UNCERTAINTIES CALCULATOR ###

# @producer(
#     uses={
#         'event',
#     },
#     produces={
#         f"PuppiMET_recoil_unc_{coord}_{shift}"
#         for coord in ("pt", "phi")
#         for shift in ("RespUp", "RespDown", "ResolUp", "ResolDown")
#     },
#     mc_only=True,
# )
# def DY_pTll_recoil_unc(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
#     shifts = []
#     shifts=[*shifts,"RespUp", "RespDown", "ResolUp", "ResolDown"] 
    
#     gen = events.GenPart

#     zeros_array =  ak.zeros_like(gen)
    
#     pdgId = abs(gen.pdgId)
#     sFlag = gen.statusFlags
#     # Pick up every charged lepton and neutrino which is “fromHardProcess” and stable
#     first_condition = ((pdgId >= 11) & (pdgId <= 16) & (sFlag >> 8 & 1) & (gen.status == 1))
#     #Also pick up anything with bit 10: “isDirectHardProcessTauDecayProduct”, these are pions from hadronic Tau decays
#     second_condition = (sFlag >> 10 & 1)
#     GenObj_mask = (first_condition | second_condition)
    
#     pdgId_neutrinos = ((pdgId == 12) | (pdgId == 14) | (pdgId == 16))
#     VisObj_mask = (GenObj_mask & (~pdgId_neutrinos))
    
#     gen_obj = ak.where(GenObj_mask, gen, zeros_array)
#     vis_obj = ak.where(VisObj_mask, gen, zeros_array)
    
#     GenObj = get_lep_p4(gen_obj)

#     VisObj = get_lep_p4(vis_obj)
    
    
#     def ptetaphim_to_xyzt(pt, eta, phi, mass):
#         px = pt * np.cos(phi)
#         py = pt * np.sin(phi)
#         pz = pt * np.sinh(eta)
#         # Energy from E^2 = (p*c)^2 + (m*c^2)^2, here c=1
#         E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
#         return px, py, pz, E

#     def xyzt_to_ptetaphim(px, py, pz, E):
#         pt = np.sqrt(px**2 + py**2)
#         # avoid divide-by-zero
#         eta = np.arcsinh(pz / pt)
#         phi = np.arctan2(py, px)
#         mass = np.sqrt(np.abs(E**2 - (px**2 + py**2 + pz**2)))
#         return pt, eta, phi, mass

#     # Suppose GenObj is a nested PtEtaPhiMLorentzVectorArray of shape
#     # [num_events] -> variable-length list of 4-vectors.

#     px, py, pz, E = ptetaphim_to_xyzt(GenObj.pt, GenObj.eta, GenObj.phi, GenObj.mass)

#     # Sum within each event (axis=1). This collapses each sublist into one vector.
#     px_sum  = ak.sum(px, axis=1)
#     py_sum  = ak.sum(py, axis=1)
#     pz_sum  = ak.sum(pz, axis=1)
#     E_sum   = ak.sum(E, axis=1)

#     # Convert back to PtEtaPhiM
#     pt_sum, eta_sum, phi_sum, mass_sum = xyzt_to_ptetaphim(px_sum, py_sum, pz_sum, E_sum)

#     # Finally zip it up again as a PtEtaPhiMLorentzVector record (if desired)
#     GenObj_sum = ak.zip(
#         {
#             "pt": pt_sum,
#             "eta": eta_sum,
#             "phi": phi_sum,
#             "mass": mass_sum,
#         },
#         with_name="PtEtaPhiMLorentzVector",
#     )
#     # Suppose VisObj is a nested PtEtaPhiMLorentzVectorArray of shape
#     # [num_events] -> variable-length list of 4-vectors.

#     px, py, pz, E = ptetaphim_to_xyzt(VisObj.pt, VisObj.eta, VisObj.phi, VisObj.mass)

#     # Sum within each event (axis=1). This collapses each sublist into one vector.
#     px_sum  = ak.sum(px, axis=1)
#     py_sum  = ak.sum(py, axis=1)
#     pz_sum  = ak.sum(pz, axis=1)
#     E_sum   = ak.sum(E, axis=1)

#     # Convert back to PtEtaPhiM
#     pt_sum, eta_sum, phi_sum, mass_sum = xyzt_to_ptetaphim(px_sum, py_sum, pz_sum, E_sum)

#     # Finally zip it up again as a PtEtaPhiMLorentzVector record (if desired)
#     VisObj_sum = ak.zip(
#         {
#             "pt": pt_sum,
#             "eta": eta_sum,
#             "phi": phi_sum,
#             "mass": mass_sum,
#         },
#         with_name="PtEtaPhiMLorentzVector",
#     )
#     # Adding necessary empty eta and mass fields to create PuppiMET Lorentz vector 
#     events["PuppiMET"] = ak.with_field(
#         events["PuppiMET"],
#         ak.zeros_like(events["PuppiMET"].pt),
#         "eta"
#     )
#     events["PuppiMET"] = ak.with_field(
#         events["PuppiMET"],
#         ak.zeros_like(events["PuppiMET"].pt),
#         "mass"
#     )
#     METvector = get_lep_p4(events.PuppiMET)
    

#     H = -METvector - VisObj_sum
#     dPhi_Z_H = H.phi - GenObj_sum.phi
#     H_para = H.pt * np.cos(dPhi_Z_H)
#     H_perp = H.pt * np.sin(dPhi_Z_H)
#     ch_str = self.config_inst.channels.names()[0]
#     hcand = events[f'hcand_{ch_str}']

#     order      = "LO"
#     njet       = events.n_jets 
#     ptll       = flat_np_view(hcand.pt)
#     var_para   = "Hpara"
#     var_perp   = "Hperp"
#     val_para   = H_para
#     val_perp   = H_perp
#     #Prepare a tuple with the inputs of the correction evaluator
#     recoil_Hpara_args = lambda syst : (order,njet,ptll,var_para,val_para,syst)

#     recoil_Hperp_args = lambda syst : (order,njet,ptll,var_perp,val_perp,syst)
    
#     #Loop over the shifts and calculate for each shift muon scale factor
#     DY_pTll_recoil_unc_pt_values = {}
#     DY_pTll_recoil_unc_phi_values = {}
#     for the_shift in shifts:
#         H_para_new = self.UncertaintyCorr.evaluate(*recoil_Hpara_args(the_shift))
#         H_perp_new = self.UncertaintyCorr.evaluate(*recoil_Hperp_args(the_shift))

#         H_pt_new = np.sqrt(H_para_new**2 + H_perp_new**2)
#         dPhi_Z_H_new = np.atan2(H_perp_new, H_para_new) 
#         H_phi_new =  dPhi_Z_H_new + GenObj_sum.phi
#         H_more_then_pi = H_phi_new > np.pi
#         H_less_then_minus_pi = H_phi_new < -np.pi
        
#         H_phi_new = ak.where(H_more_then_pi, H_phi_new - 2*np.pi, H_phi_new)
#         H_phi_new = ak.where(H_less_then_minus_pi, H_phi_new + 2*np.pi, H_phi_new)
        
#         H_new = ak.zip(
#             {
#                 "pt": H_pt_new,
#                 "eta": METvector.eta,
#                 "phi": H_phi_new,
#                 "mass": METvector.mass,
#             },
#             with_name="PtEtaPhiMLorentzVector",
#         )
#         METvector_new = -H_new - VisObj_sum  

#         DY_pTll_recoil_unc_pt_values[the_shift] = METvector_new.pt
#         DY_pTll_recoil_unc_phi_values[the_shift] = METvector_new.phi
        
#     for the_shift in shifts: 
#         events = set_ak_column_f32(events, f"PuppiMET_recoil_unc_pt_{the_shift}", DY_pTll_recoil_unc_pt_values[the_shift])
#         events = set_ak_column_f32(events, f"PuppiMET_recoil_unc_phi_{the_shift}", DY_pTll_recoil_unc_phi_values[the_shift])

#     return events

# @DY_pTll_recoil_unc.requires
# def DY_pTll_recoil_unc_requires(self: Producer, reqs: dict) -> None:
#     if "external_files" in reqs:
#         return
    
#     from columnflow.tasks.external import BundleExternalFiles
#     reqs["external_files"] = BundleExternalFiles.req(self.task)

# @DY_pTll_recoil_unc.setup
# def DY_pTll_recoil_unc_setup(
#     self: Producer,
#     reqs: dict,
#     inputs: dict,
#     reader_targets: InsertableDict,
# ) -> None:
#     bundle = reqs["external_files"]
#     import correctionlib
#     correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
#     correction_set = correctionlib.CorrectionSet.from_string(
#         bundle.files.DY_pTll_recoil_corr.load(formatter="gzip").decode("utf-8"),
#     )
#     self.UncertaintyCorr = correction_set["Recoil_correction_Uncertainty"]