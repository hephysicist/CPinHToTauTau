# coding: utf-8

"""
Producers for the TauTheDifference BDT for signal vs. background separation of the Higgs CP analysis.
"""
import law
import os
import sys 

sys.path.append(os.path.dirname(os.path.realpath(__file__))) #This is needed to import fastmtt_cpp.so
import fastmtt_cpp

import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, dev_sandbox
from columnflow.util import DotDict, InsertableDict
from columnflow.columnar_util import EMPTY_FLOAT, set_ak_column, flat_np_view

np = maybe_import("numpy")
ak = maybe_import("awkward")
pd = maybe_import("pandas")
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)
set_ak_column_i32 = functools.partial(set_ak_column, value_type=np.int32)


@producer(
    uses={
        "event",
        "hcand*", "PuppiMET*",
    },
    produces = {"hcand*"}
)
def fast_mtt(
    self: Producer,
    events: ak.Array,
    **kwargs,
) -> ak.Array:
    """
    Returns the  the result of fast_mtt algorithm
    """

    ch_str = self.config_inst.channels.names()[0]
    ch_obj = self.config_inst.x.ch_objects[ch_str]
    hcand = events[f'hcand_{ch_str}']
    met = events.PuppiMET
    
    print('Running fastMTT (be patient)')

    # steering parameters
    verbosity = False
    delta = 1.0/1.15 # regularization parameter delta
    reg_order = 6.0  # regularization parameter order
    mX = 125.10 # Higgs mass
    widthX = 2.5 # window
    # if args.sample=='dy':
    # mX = 91.2 # Z boson mass
    # widthX = 4.0 # window

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
    
    
    fast_mtt_res = fastmtt_cpp.fastmtt_cpp(*features_list)
    fast_mtt_unflattened = {}
    for the_var, the_arr in fast_mtt_res.items():
        fast_mtt_unflattened[the_var] = ak.unflatten(the_arr, ak.num(hcand.mass,axis=1)) 
    hcand['fastMTT'] = ak.zip(fast_mtt_unflattened)
    events = set_ak_column(events, f'hcand_{ch_str}', hcand)
    return events


