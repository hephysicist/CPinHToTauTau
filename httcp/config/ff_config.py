# coding: utf-8

"""
Definition of categories.
"""

import order as od

from law.util import DotDict
from columnflow.util import maybe_import
np = maybe_import("numpy")
import itertools
from columnflow.config_util import add_shift_aliases


def add_ff_config(config: od.Config,
                   channel = None) -> None:
    

    config.x.fake_factor_method = DotDict.wrap({
        "axes": {'tau_pt': {
                    'var_route' : [f'hcand_{channel}', 'lep1', 'pt'],
                    'ax_str'    : 'Variable([20,25,30,35,40,50,60,80,300], name="tau_pt", label="Tau pt", underflow=False, overflow=False)',
                    },
                'tau_dm_pnet': {
                    'var_route' : [f'hcand_{channel}', 'lep1', 'decayModePNet'],
                    'ax_str'    : 'IntCategory([0,1,2,10,11], name="tau_dm_pnet", label="Tau PNet decayMode")',
                },
                "n_jets": {
                    'var_route' : ['n_jets'],
                    'ax_str'    : 'Integer(0, 3, name="n_j", label="Number of jets",underflow=False, overflow=False)',
                },
        },
        "columns" : ['ff_wj','ff_qcd'],
        "shifts"  : {
            'tau_dm_pnet'   : [0,1,2,10],
            "n_j"           : [0,1,2],        
            "shift_name"    : ["up", "down"],
            },
        })
    shifts = config.x.fake_factor_method.shifts
    ff_type = config.x.fake_factor_method.columns
    
    shift_tuple = itertools.product(ff_type,
                                    shifts.tau_dm_pnet,
                                    shifts.n_j)
    
    syst_names = []
    for i, (the_name, the_dm, the_nj) in enumerate(shift_tuple):
        syst_name = '_'.join(('dm',str(the_dm),
                        'nj',str(the_nj)))
        syst_names.append('_'.join((the_name, syst_name)))
        
        config.add_shift(name='_'.join((the_name, syst_name, "up")) , id=1000 + 2 * i, type="shape", tags={"ff"})
        config.add_shift(name= '_'.join((the_name, syst_name, "down")), id=1001 + 2 * i, type="shape", tags={"ff"})
        add_shift_aliases(config, 
                          '_'.join((the_name, syst_name)),
                          {
                              f"ff_weight.{the_name}"    : "ff_weight.{name}",
                          },
        )
    config.x.ff_syst_names = syst_names