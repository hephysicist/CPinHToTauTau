# coding: utf-8

"""
A producer to apply FastMTT to the lepton candidates in the hcand columns.
This reconstructs the full 4-momenta of the tau leptons via a scan over the 
di-tau mass hypothesis, using the MET as a constraint on the neutrino system.

Use cases :
-> corrected tau pTs -> phicp
-> corrected di-Tau mass -> TauTheDifference BDT for signal vs. background separation of the Higgs CP analysis
"""

import functools
import law
import os
import sys
import order as od
from typing import Optional

sys.path.append(os.path.dirname(os.path.realpath(__file__))) #This is needed to import fastmtt_cpp.so
import fastmtt_cpp

from columnflow.production import Producer, producer
from columnflow.production.util import attach_coffea_behavior

from columnflow.util import maybe_import, dev_sandbox
from columnflow.util import DotDict, InsertableDict
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, flat_np_view
from columnflow.columnar_util import optional_column as optional

from httcp.util import get_lep_p4

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
pd = maybe_import("pandas")
cnmn = maybe_import("coffea.nanoevents.methods.nanoaod")

# helpers
set_ak_column_f32 = functools.partial(set_ak_column, value_dc_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_dc_type=np.int32)

logger = law.logger.get_logger(__name__)

@producer(
    uses={
        # nano columns
        'event',
        'hcand_*',
        'PuppiMET*', #here : change from PuppiMET.pt, PuppiMET.phi, PuppiMET.covXX, PuppiMET.covXY, PuppiMET.covYY
    },
    produces={
        # new columns
        'hcand_*',
    },
)

def fastMTT(
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

    ch_str = self.config_inst.channels.names()[0] # Gives list of strings with current channel ID e.g 'mutau'
    ch_obj = self.config_inst.x.ch_objects[ch_str]

    hcand = events[f'hcand_{ch_str}']
    met = events.PuppiMET

    print('Running fastMTT')


    # Steering Parameters for the FastMTT algorithm
    verbosity = False # True = prints fastMTT infos in terminal
    delta = 1.0/1.15 # regularization parameter delta
    reg_order = 6.0  # regularization parameter order
    mX = 125.10 # Higgs mass
    widthX = 2.5 # window
    #if args.sample=='dy':
    #mX = 91.2 # Z boson mass
    #widthX = 4.0 # window

    # Prepare input features
    N = len(flat_np_view(hcand.lep0.pt, axis=1))
    lep_features = ['pt','eta','phi','mass']
    
    features_list = [int(N)]
    
    for lep_str in ['lep0', 'lep1']:
        lep = getattr(hcand, lep_str)
        for the_var in lep_features:
            features_list.append(flat_np_view(getattr(lep, the_var), axis=1).astype(np.float64))
            ####################################
            # decay_type = 0 : tau -> electron #
            #            = 1 : tau -> muon     #
            #            = 2 : tau -> hadrons  #
            ####################################
        decay_type = np.ones_like(flat_np_view(lep.pt, axis=1),dtype=np.int32)
        if ch_obj[lep_str] == 'Electron'    : decay_type = np.right_shift(1, decay_type)
        if ch_obj[lep_str] == 'Tau'         : decay_type = np.left_shift(1, decay_type)
        features_list.append(decay_type)
            
    met_features = ['pt','phi','covXX','covXY','covYY']
    for the_var in met_features:
        features_list.append(flat_np_view(getattr(met, the_var)).astype(np.float64))
    
    fastmtt_kwargs = [verbosity,delta,reg_order,mX,widthX]
    for the_var in fastmtt_kwargs:
        features_list.append(the_var)
    
    
    # === Run FastMTT ===
    fast_mtt_res = fastmtt_cpp.fastmtt_cpp(*features_list)
    fast_mtt_unflattened = {}
    for the_var, the_arr in fast_mtt_res.items():
        fast_mtt_unflattened[the_var] = ak.unflatten(the_arr, ak.num(hcand.mass,axis=1)) 
    MTT_cpp_rest = ak.zip(fast_mtt_unflattened) # -> mass, x1, x1_BW, x1_cons, x2, x2_BW, x2_cons

    ''' ATTENTION : x1_BW and x2_BW denote Boson Mass Window, not Breit Wigner ! 
    
        (Breit Wigner can be used in apply_fastMTT.py and 
        ~/CPinHToTauTau/modules/ClassicsSVfit/python/FastMTT.py l.37)'''


    # Extract reconstructed MTT Higgs mass and decay products properties
    hcand_features = ['pt', 'eta', 'phi'] #no 'px', 'py', 'pz', in lep so far (put can be used in apply_fastMTT.py)
    leptons_corrected = {}
    
    # Replace possible NaNs in lep0 mass with SM values
    lep0_mass_defaults = {
        'mutau': 0.10566,   # muon mass in GeV
        'etau': 0.000511,   # electron mass in GeV
        'tautau': 1.77686   # tau mass in GeV 
        }
    lep0_mass_threshold = 1e-7

    for lep_str in ['lep0', 'lep1']:
        lep = getattr(hcand, lep_str)

        modifier_key = 'x1' if lep_str == 'lep0' else 'x2'
        modifier = MTT_cpp_rest[modifier_key]

        lep_dict = {}

        for hcand_feature in hcand_features:
            lep_dict[hcand_feature] = getattr(lep, hcand_feature) / modifier

        # Handle mass separately, as it is not a Lorentz vector component
        lep_p4 = get_lep_p4(lep) # visible ['pt', 'eta', 'phi', 'mass']
        lep_energy = lep_p4.energy
        lep_dict['energy'] = lep_energy / modifier

        lep_p4_mtt = ak.zip(lep_dict,
                    with_name="PtEtaPhiELorentzVector",
                    behavior=coffea.nanoevents.methods.vector.behavior)

        if lep_str == 'lep0':
            lep_dict['mass'] = ak.where(lep_p4_mtt.mass < lep0_mass_threshold,
                                            lep0_mass_defaults.get(ch_str, 0.0),
                                            lep_p4_mtt.mass)
        else:
            lep_dict['mass'] = lep_p4_mtt.mass
            
        # Prepare the reconstructed tau leptons into new hcand columns | part 1/2
        leptons_corrected[lep_str] = ak.zip(lep_dict)

    # Prepare the reconstructed tau leptons into new hcand columns | part 2/2
    fastMTT = ak.zip({
            'lep0': leptons_corrected['lep0'],
            'lep1': leptons_corrected['lep1'],
            'mass': fast_mtt_unflattened['mass'],
        })

    # Attach the new fastMTT entry to the corresponding hcand column
    hcand['fastMTT'] = fastMTT
    events = set_ak_column(events, f'hcand_{ch_str}', hcand)

    return events
