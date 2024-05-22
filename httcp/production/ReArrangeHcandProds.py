import os

from typing import Optional
from columnflow.util import maybe_import
from columnflow.production import Producer, producer
#from httcp.production.ReconstructPi0 import reconstructPi0

from httcp.util import getGenTauDecayMode

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")



def convert_to_coffea_p4(zipped_item):
    return ak.zip(
        zipped_item,
        with_name = "PtEtaPhiMLorentzVector",
        behavior  = coffea.nanoevents.methods.vector.behavior,
    )

def getMaxEtaTauStrip(pt):
    temp = 0.20 * np.power(pt, -0.66)
    ref1 = ak.where(temp > 0.15, 0.15, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2

def getMaxPhiTauStrip(pt):
    temp = 0.35 * np.power(pt, -0.71)
    ref1 = ak.where(temp > 0.30, 0.30, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2


def reconstructPi0(hcandp4, photons, method: Optional[str] = "simpleIC"):
    photons = ak.with_name(photons, "PtEtaPhiMLorentzVector")
    photons_sorted_pt_indices = ak.argsort(photons.pt, ascending=False)
    photons = photons[photons_sorted_pt_indices]

    p4_pi0 = None

    if method == "simpleIC":
        photons_px = ak.sum(photons.px, axis=1)
        photons_py = ak.sum(photons.py, axis=1)
        photons_pt = np.sqrt(photons_px ** 2 + photons_py ** 2)

        #pt_pi0    = photons[:, 0:1].pt
        pt_pi0    = photons_pt
        eta_pi0   = photons[:, 0:1].eta
        phi_pi0   = photons[:, 0:1].phi
        pdgid_pi0 = ak.values_astype(111 * ak.ones_like(eta_pi0), "int64")
        mass_pi0  = 0.135 * ak.ones_like(eta_pi0)
        
        p4_pi0 = convert_to_coffea_p4({
            "pt"    : pt_pi0,
            "eta"   : eta_pi0,
            "phi"   : phi_pi0,
            "mass"  : mass_pi0,
            "pdgId" : pdgid_pi0,
        })
        
    elif method == "simpleMB":
        has_atleast_one_photon = ak.num(photons.pt, axis=1) > 0
        hcandp4 = ak.where(has_atleast_one_photon, hcandp4, hcandp4[:,:0])
        photons_p4 = ak.where(has_atleast_one_photon,
                              photons_p4[:,0:1],
                              photons_p4[:,:0])
        
        deta_photons_hcand = (photons_p4).metric_table(hcandp4, metric = lambda a,b: np.abs(a.eta - b.eta))
        dphi_photons_hcand = (photons_p4).metric_table(hcandp4, metric = lambda a,b: np.abs(a.delta_phi(b)))
        
        maxeta_photons = getMaxEtaTauStrip(photons_p4.pt)
        maxphi_photons = getMaxPhiTauStrip(photons_p4.pt)

        mask_photons = ((np.abs(deta_photons_hcand) < maxeta_photons)
                        & (np.abs(dphi_photons_hcand) < maxphi_photons))
        
        strip_photons_p4 = convert_to_coffea_p4(
            {
                "pt"   : photons_p4.pt[mask_photons],
                "eta"  : photons_p4.eta[mask_photons],
                "phi"  : photons_p4.phi[mask_photons],
                "mass" : photons_p4.mass[mask_photons],
            }
        )

        has_one_photon = ak.num(strip_photons_p4.pt, axis=1) == 1
        
        strip_photons_p4_pair = ak.combinations(strip_photons_p4, 2, axis=1)
        strip_photons_p4_pair_0, strip_photons_p4_pair_1 = ak.unzip(strip_photons_p4_pair)
        strip_photons_mass = (strip_photons_p4_pair_0 + strip_photons_p4_pair_1).mass
        mass_val = np.abs(strip_photons_mass - pi0RecoM)
        mass_sorted_idx = ak.argsort(mass_val, axis=1)
        strip_photons_p4_pair_sorted = strip_photons_p4_pair[mass_sorted_idx]
        mass_mask = mass_val < 2 * pi0RecoW
        evt_mask_no_pair = ak.sum(mass_mask, axis=1) == 0
        strip_photons_p4_pair_sorted_pass_mass = strip_photons_p4_pair_sorted[mass_mask]
        strip_photons_p4_mass_selected = ak.concatenate([strip_photons_p4_pair_sorted_pass_mass["0"][:,:1],
                                                         strip_photons_p4_pair_sorted_pass_mass["1"][:,:1]], axis=1)

        sel_strip_photons_p4 = ak.where(has_one_photon,
                                        strip_photons_p4,
                                        ak.where(evt_mask_no_pair,
                                                 strip_photons_p4[:,:1],
                                                 strip_photons_p4_mass_selected)
                                    )
        
        sel_strip_pizero_p4 = ak.from_regular(ak.sum(sel_strip_photons_p4, axis=-1)[:,None])

        p4_pi0 = sel_strip_pizero_p4


    return p4_pi0


def getpions(decay_gentau: ak.Array) -> ak.Array :
    ispion_pos = lambda prod: ((prod.pdgId ==  211) | (prod.pdgId ==  321))
    ispion_neg = lambda prod: ((prod.pdgId == -211) | (prod.pdgId == -321))

    pions_tau  = decay_gentau[(ispion_pos(decay_gentau) | ispion_neg(decay_gentau))]

    has_three_pions = lambda prod : ak.sum((ispion_pos(prod) | ispion_neg(prod)), axis=1) == 3
    has_two_pos_one_neg_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 2) & (ak.sum(ispion_neg(prod), axis=1) == 1)
    has_two_neg_one_pos_pions = lambda prod : has_three_pions(prod) & (ak.sum(ispion_pos(prod), axis=1) == 1) & (ak.sum(ispion_neg(prod), axis=1) == 2)
    get_sorted_pion_indices   = lambda pions: ak.where(has_two_pos_one_neg_pions(pions),
                                                       ak.argsort(pions.pdgId, ascending=True),
                                                       ak.where(has_two_neg_one_pos_pions(pions),
                                                                ak.argsort(pions.pdgId, ascending=False),
                                                                ak.local_index(pions.pdgId)))
    
    pions_tau_sorted_indices = get_sorted_pion_indices(pions_tau)
    sorted_pions_tau = pions_tau[pions_tau_sorted_indices]

    return sorted_pions_tau


def getphotons(decay_tau: ak.Array) -> ak.Array :
    isphoton   = lambda prod: (prod.pdgId == 22)
    photons = decay_tau[isphoton(decay_tau)]
    return photons


def getgenpizeros(decay_gentau: ak.Array) -> ak.Array :
    #ispizero = lambda col: (np.abs(col.pdgId) == 111) | (np.abs(col.pdgId) == 311) | (np.abs(col.pdgId) == 130) | (np.abs(col.pdgId) == 310)
    ispizero = lambda col: ((col.pdgId == 111) 
                            | (col.pdgId == 311) 
                            | (col.pdgId == 130) 
                            | (col.pdgId == 310))
    pizeros_tau = decay_gentau[ispizero(decay_gentau)]

    return pizeros_tau


@producer(
    uses={
        "channel_id", 
        "hcand.pt", "hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode",
        "hcandprod.pt", "hcandprod.eta", "hcandprod.phi", "hcandprod.mass","hcandprod.pdgId",
    },
)
def reArrangeDecayProducts(
        self: Producer,
        events: ak.Array,
        **kwargs
) :
    hcand      = events.hcand
    hcandprod  = events.hcandprod

    hcand1     = hcand[:, 0:1]
    hcand2     = hcand[:, 1:2]
    hcand1prod = ak.firsts(hcandprod[:,0:1], axis=1)
    hcand2prod = ak.firsts(hcandprod[:,1:2], axis=1)

    hcand1prod_photons = getphotons(hcand1prod)
    hcand2prod_photons = getphotons(hcand2prod)
    
    hcand1prod_pions = getpions(hcand1prod)
    hcand2prod_pions = getpions(hcand2prod)
    
    # hcand1 and its decay products
    p4_hcand1     = ak.with_name(hcand1, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi  = ak.with_name(hcand1prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi0 = reconstructPi0(p4_hcand1, hcand1prod_photons)

    # hcand2 and its decay products
    p4_hcand2     = ak.with_name(hcand2, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi  = ak.with_name(hcand2prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi0 = reconstructPi0(p4_hcand2, hcand2prod_photons)

    hcand1AndProds = ak.concatenate([p4_hcand1, p4_hcand1_pi, p4_hcand1_pi0], axis=1)
    hcand2AndProds = ak.concatenate([p4_hcand2, p4_hcand2_pi, p4_hcand2_pi0], axis=1)

    #from IPython import embed; embed()

    return events, {"p4h1"       : p4_hcand1, 
                    "p4h1pi"     : p4_hcand1_pi, 
                    "p4h1pi0"    : p4_hcand1_pi0, 
                    "p4h2"       : p4_hcand2, 
                    "p4h2pi"     : p4_hcand2_pi, 
                    "p4h2pi0"    : p4_hcand2_pi0}

@producer(
    uses={
        "channel_id",
        "GenTau.*",
        "GenTauProd.*"
    },
    #produces={
    #    "GenTau.decayMode",
    #},
    mc_only=True,
)
def reArrangeGenDecayProducts(
        self: Producer,
        events: ak.Array,
        **kwargs
) -> tuple[ak.Array, dict] :
    ghcand       = events.GenTau
    ghcandprod   = events.GenTauProd
    #from IPython import embed; embed()
    #ghcandprod_indices = ghcand.distinctChildrenIdxG
    #ghcandprod         = events.GenPart._apply_global_index(ghcandprod_indices)

    #hcandprod_dm = getGenTauDecayMode(ghcandprod)
    #events = set_ak_column(events, "GenTau.decayMode", hcandprod_dm)

    hcand1     = ghcand[:, 0:1]
    hcand2     = ghcand[:, 1:2]
    hcand1prod = ak.firsts(ghcandprod[:,0:1], axis=1)
    hcand2prod = ak.firsts(ghcandprod[:,1:2], axis=1)

    hcand1prod_pions = getpions(hcand1prod)
    hcand2prod_pions = getpions(hcand2prod)

    hcand1prod_pizeros = getgenpizeros(hcand1prod)
    hcand2prod_pizeros = getgenpizeros(hcand2prod)
        
    # hcand1 and its decay products
    p4_hcand1     = ak.with_name(hcand1, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi  = ak.with_name(hcand1prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand1_pi0 = ak.with_name(hcand1prod_pizeros, "PtEtaPhiMLorentzVector")

    # hcand2 and its decay products
    p4_hcand2     = ak.with_name(hcand2, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi  = ak.with_name(hcand2prod_pions, "PtEtaPhiMLorentzVector")
    p4_hcand2_pi0 = ak.with_name(hcand2prod_pizeros, "PtEtaPhiMLorentzVector")

    hcand1AndProds = ak.concatenate([p4_hcand1, p4_hcand1_pi, p4_hcand1_pi0], axis=1)
    hcand2AndProds = ak.concatenate([p4_hcand2, p4_hcand2_pi, p4_hcand2_pi0], axis=1)

    #from IPython import embed; embed()

    return events, {"p4h1"        : p4_hcand1, 
                    "p4h1pi"      : p4_hcand1_pi, 
                    "p4h1pi0"     : p4_hcand1_pi0, 
                    "p4h2"        : p4_hcand2, 
                    "p4h2pi"      : p4_hcand2_pi, 
                    "p4h2pi0"     : p4_hcand2_pi0}
