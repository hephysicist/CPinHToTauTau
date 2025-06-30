import functools

from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior
from columnflow.selection.util import sorted_indices_from_mask

from MSSM_H_tt.util import get_lep_p4, get_p2

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
cl = maybe_import("correctionlib")
warn = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


# "name": "Recoil_correction_Uncertainty",
#    "description": "Various values needed to recompute the PuppiMET of single-boson processes like DY. The values were derived from a mu-mu control region, while validation and uncertainties are obtained from e-e control region. 'Recoil_correction_Uncertainty' gives uncertainties which are to be applies in addition to the nominal corrections. Note that a different input is required. The uncertainties were derived based on corrections applied with the QuantileMapHist method.",
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
#      "description": "The variable name you are giving the input for: 'Hpara', 'Hperp' (string). Note that this is a DIFFERENT vairable from Upara/Uperp! The output will be for the same kind of variable."
#     },
#     {
#      "name": "val",
#      "type": "real",
#      "description": "Input value of either Hpara or Hperp"
#     },
#     {
#      "name": "syst",
#      "type": "string",
#      "description": "Systematic variations on Response and Resolution: 'RespUp', 'RespDown', 'ResolUp', 'ResolDown'"
#     }


### DY PTLL RECOIL UNCERTAINTIES CALCULATOR ###

@producer(
    uses={
        'event',
    },
    produces={
        f"PuppiMET_recoil_unc_{coord}_{shift}"
        for coord in ("pt", "phi")
        for shift in ("RespUp", "RespDown", "ResolUp", "ResolDown")
    },
    mc_only=True,
)
def DY_pTll_recoil_unc(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    
    shifts = []
    shifts=[*shifts,"RespUp", "RespDown", "ResolUp", "ResolDown"] 
    
    gen = events.GenPart

    zeros_array =  ak.zeros_like(gen)
    
    pdgId = abs(gen.pdgId)
    sFlag = gen.statusFlags
    # Pick up every charged lepton and neutrino which is “fromHardProcess” and stable
    first_condition = ((pdgId >= 11) & (pdgId <= 16) & (sFlag >> 8 & 1) & (gen.status == 1))
    #Also pick up anything with bit 10: “isDirectHardProcessTauDecayProduct”, these are pions from hadronic Tau decays
    second_condition = (sFlag >> 10 & 1)
    GenObj_mask = (first_condition | second_condition)
    
    pdgId_neutrinos = ((pdgId == 12) | (pdgId == 14) | (pdgId == 16))
    VisObj_mask = (GenObj_mask & (~pdgId_neutrinos))
    
    gen_obj = ak.where(GenObj_mask, gen, zeros_array)
    vis_obj = ak.where(VisObj_mask, gen, zeros_array)
    
    GenObj = get_lep_p4(gen_obj)

    VisObj = get_lep_p4(vis_obj)
    
    
    def ptetaphim_to_xyzt(pt, eta, phi, mass):
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        # Energy from E^2 = (p*c)^2 + (m*c^2)^2, here c=1
        E = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
        return px, py, pz, E

    def xyzt_to_ptetaphim(px, py, pz, E):
        pt = np.sqrt(px**2 + py**2)
        # avoid divide-by-zero
        eta = np.arcsinh(pz / pt)
        phi = np.arctan2(py, px)
        mass = np.sqrt(np.abs(E**2 - (px**2 + py**2 + pz**2)))
        return pt, eta, phi, mass

    # Suppose GenObj is a nested PtEtaPhiMLorentzVectorArray of shape
    # [num_events] -> variable-length list of 4-vectors.

    px, py, pz, E = ptetaphim_to_xyzt(GenObj.pt, GenObj.eta, GenObj.phi, GenObj.mass)

    # Sum within each event (axis=1). This collapses each sublist into one vector.
    px_sum  = ak.sum(px, axis=1)
    py_sum  = ak.sum(py, axis=1)
    pz_sum  = ak.sum(pz, axis=1)
    E_sum   = ak.sum(E, axis=1)

    # Convert back to PtEtaPhiM
    pt_sum, eta_sum, phi_sum, mass_sum = xyzt_to_ptetaphim(px_sum, py_sum, pz_sum, E_sum)

    # Finally zip it up again as a PtEtaPhiMLorentzVector record (if desired)
    GenObj_sum = ak.zip(
        {
            "pt": pt_sum,
            "eta": eta_sum,
            "phi": phi_sum,
            "mass": mass_sum,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    # Suppose VisObj is a nested PtEtaPhiMLorentzVectorArray of shape
    # [num_events] -> variable-length list of 4-vectors.

    px, py, pz, E = ptetaphim_to_xyzt(VisObj.pt, VisObj.eta, VisObj.phi, VisObj.mass)

    # Sum within each event (axis=1). This collapses each sublist into one vector.
    px_sum  = ak.sum(px, axis=1)
    py_sum  = ak.sum(py, axis=1)
    pz_sum  = ak.sum(pz, axis=1)
    E_sum   = ak.sum(E, axis=1)

    # Convert back to PtEtaPhiM
    pt_sum, eta_sum, phi_sum, mass_sum = xyzt_to_ptetaphim(px_sum, py_sum, pz_sum, E_sum)

    # Finally zip it up again as a PtEtaPhiMLorentzVector record (if desired)
    VisObj_sum = ak.zip(
        {
            "pt": pt_sum,
            "eta": eta_sum,
            "phi": phi_sum,
            "mass": mass_sum,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    # Adding necessary empty eta and mass fields to create PuppiMET Lorentz vector 
    events["PuppiMET"] = ak.with_field(
        events["PuppiMET"],
        ak.zeros_like(events["PuppiMET"].pt),
        "eta"
    )
    events["PuppiMET"] = ak.with_field(
        events["PuppiMET"],
        ak.zeros_like(events["PuppiMET"].pt),
        "mass"
    )
    METvector = get_lep_p4(events.PuppiMET)
    

    H = -METvector - VisObj_sum
    dPhi_Z_H = H.phi - GenObj_sum.phi
    H_para = H.pt * np.cos(dPhi_Z_H)
    H_perp = H.pt * np.sin(dPhi_Z_H)
    ch_str = self.config_inst.channels.names()[0]
    hcand = events[f'hcand_{ch_str}']

    order      = "LO"
    njet       = events.n_jets 
    ptll       = flat_np_view(hcand.pt)
    var_para   = "Hpara"
    var_perp   = "Hperp"
    val_para   = H_para
    val_perp   = H_perp
    #Prepare a tuple with the inputs of the correction evaluator
    recoil_Hpara_args = lambda syst : (order,njet,ptll,var_para,val_para,syst)

    recoil_Hperp_args = lambda syst : (order,njet,ptll,var_perp,val_perp,syst)
    
    #Loop over the shifts and calculate for each shift muon scale factor
    DY_pTll_recoil_unc_pt_values = {}
    DY_pTll_recoil_unc_phi_values = {}
    for the_shift in shifts:
        H_para_new = self.UncertaintyCorr.evaluate(*recoil_Hpara_args(the_shift))
        H_perp_new = self.UncertaintyCorr.evaluate(*recoil_Hperp_args(the_shift))

        H_pt_new = np.sqrt(H_para_new**2 + H_perp_new**2)
        dPhi_Z_H_new = np.atan2(H_perp_new, H_para_new) 
        H_phi_new =  dPhi_Z_H_new + GenObj_sum.phi
        H_more_then_pi = H_phi_new > np.pi
        H_less_then_minus_pi = H_phi_new < -np.pi
        
        H_phi_new = ak.where(H_more_then_pi, H_phi_new - 2*np.pi, H_phi_new)
        H_phi_new = ak.where(H_less_then_minus_pi, H_phi_new + 2*np.pi, H_phi_new)
        
        H_new = ak.zip(
            {
                "pt": H_pt_new,
                "eta": METvector.eta,
                "phi": H_phi_new,
                "mass": METvector.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        METvector_new = -H_new - VisObj_sum  

        DY_pTll_recoil_unc_pt_values[the_shift] = METvector_new.pt
        DY_pTll_recoil_unc_phi_values[the_shift] = METvector_new.phi
        
    for the_shift in shifts: 
        events = set_ak_column_f32(events, f"PuppiMET_recoil_unc_pt_{the_shift}", DY_pTll_recoil_unc_pt_values[the_shift])
        events = set_ak_column_f32(events, f"PuppiMET_recoil_unc_phi_{the_shift}", DY_pTll_recoil_unc_phi_values[the_shift])

    return events

@DY_pTll_recoil_unc.requires
def DY_pTll_recoil_unc_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@DY_pTll_recoil_unc.setup
def DY_pTll_recoil_unc_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.DY_pTll_recoil_corr.load(formatter="gzip").decode("utf-8"),
    )
    self.UncertaintyCorr = correction_set["Recoil_correction_Uncertainty"]