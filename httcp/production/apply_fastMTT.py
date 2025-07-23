# coding: utf-8

"""
A producer to apply FastMTT to the lepton candidates in the hcand columns.
This reconstructs the full 4-momenta of the tau leptons via a scan over the 
di-tau mass hypothesis, using the MET as a constraint on the neutrino system.
"""
import functools

import law
import order as od
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column
from columnflow.columnar_util import optional_column as optional

# FastMTT implementation (precompiled C++ backend wrapped for Python)
from modules.ClassicsSVfit.python.FastMTT import FastMTT

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
cnmn = maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_dc_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_dc_type=np.int32)

logger = law.logger.get_logger(__name__)

@producer(
    uses={
        # nano columns
        'hcand_*',
        'PuppiMET.pt', 'PuppiMET.phi', 'PuppiMET.covXX', 'PuppiMET.covXY', 'PuppiMET.covYY',
    },
    produces={
        # new columns
        'hcand_*',
    },
)

def apply_fastMTT(
        self: Producer, 
        events: ak.Array,
        **kwargs
) -> ak.Array:

    def make_vector(p4):
        return ak.zip({
            'x': ak.from_regular(ak.Array(p4[:, 0:1])),
            'y': ak.from_regular(ak.Array(p4[:, 1:2])),
            'z': ak.from_regular(ak.Array(p4[:, 2:3])),
            't': ak.from_regular(ak.Array(p4[:, 3:4])),
        }, with_name='LorentzVector', behavior=coffea.nanoevents.methods.vector.behavior)


    channels = self.config_inst.channels.names() # Gives list of strings with current channel ID e.g 'mutau'

    for ch_str in channels:
        # Access the candidate collection for the current channel
        hcand_chnl = getattr(events, f'hcand_{ch_str}')

        hcand_ = ak.with_name(hcand_chnl, 'PtEtaPhiMLorentzVector')
        lep0 = ak.with_name(hcand_chnl.lep0, 'PtEtaPhiMLorentzVector')
        lep1 = ak.with_name(hcand_chnl.lep1, 'PtEtaPhiMLorentzVector')
        
        hmass_vis = (lep0 + lep1).mass

        # Map channel to FastMTT tau type
        type_map = {'mutau': 3, 'etau': 2, 'tautau': 1}
        type1_val = type_map.get(ch_str, 1)

        type1 = np.full(len(hmass_vis), type1_val, dtype=np.float32).reshape(-1, 1)
        type2 = np.full(len(hmass_vis), 1, dtype=np.float32).reshape(-1, 1)

        # Prepare the input arrays for measuredTauLeptons for FastMTT
        pt1, eta1, phi1, mass1 = map(ak.to_numpy, (lep0.pt, lep0.eta, lep0.phi, lep0.mass))
        pt2, eta2, phi2, mass2 = map(ak.to_numpy, (lep1.pt, lep1.eta, lep1.phi, lep1.mass))
        dm2 = ak.to_numpy(lep1.decayMode)
        dm1_dummy = np.zeros_like(dm2)              # dummy bc DM is constant

        metpt, metphi = map(ak.to_numpy, (events.PuppiMET.pt, events.PuppiMET.phi))
        metx, mety = metpt * np.cos(metphi), metpt * np.sin(metphi)
        
        metcov = np.stack([
            ak.to_numpy(events.PuppiMET.covXX),
            ak.to_numpy(events.PuppiMET.covXY),
            ak.to_numpy(events.PuppiMET.covXY),
            ak.to_numpy(events.PuppiMET.covYY),
        ], axis=-1).reshape(-1, 2, 2)

        # Build FastMTT input (N_events, 2 leptons, 6 features)
        measuredTauLeptons = np.stack([
            dm1_dummy, pt1, eta1, phi1, mass1, type1,
            dm2,       pt2, eta2, phi2, mass2, type2
        ], axis=1).reshape(-1, 2, 6)

        
        ## === Run FastMTT ===
        fMTT = FastMTT()
        fMTT.WhichLikelihoodPlot = -1               # -1 : likelihood plot disabled
        fMTT.CalculateUncertainties = True          # Enable uncertainty estimation

        fMTT.run(measuredTauLeptons, metx, mety, metcov)

        # Extract reconstructed MTT Higgs mass and convert to awkward array
        hmass_MTT = fMTT.mass
        hmass_MTT = ak.from_regular(ak.Array(hmass_MTT.reshape(hmass_MTT.shape[0],1)))
        
        p4_h1_reg = make_vector(fMTT.tau1P4)
        p4_h2_reg = make_vector(fMTT.tau2P4)
        
        # Replace possible NaNs in lep0 mass with SM values
        lep0_mass_defaults = {
            'mutau': 0.10566, # muon mass in GeV
            'etau': 0.000511, # electron mass in GeV
            'tautau': 1.77686 # tau mass in GeV #TODO do I need this ?
        }
        lep0_mass_4v = ak.nan_to_num(p4_h1_reg.mass, nan=lep0_mass_defaults.get(ch_str, 0.0))

        # Prepare the reconstructed tau leptons into new hcand columns
        lep0 = ak.zip(
            {
                'px': p4_h1_reg.x,
                'py': p4_h1_reg.y,
                'pz': p4_h1_reg.z,
                'pt': p4_h1_reg.pt,
                'eta': p4_h1_reg.eta,
                'phi': p4_h1_reg.phi,
                'mass': lep0_mass_4v,
            }
        )
        lep1 = ak.zip(
            {
                'px': p4_h2_reg.x,
                'py': p4_h2_reg.y,
                'pz': p4_h2_reg.z,
                'pt': p4_h2_reg.pt,
                'eta': p4_h2_reg.eta,
                'phi': p4_h2_reg.phi,
                'mass': p4_h2_reg.mass,
            }
        )

        fastMTT = ak.zip(
            {
                'lep0': lep0,
                'lep1': lep1,
                'mass': hmass_MTT,
            }
        )

        # Attach the new fastMTT entry to the corresponding hcand column
        hcand_chnl['fastMTT'] = fastMTT
        events = set_ak_column(events, f'hcand_{ch_str}', hcand_chnl)

    return events