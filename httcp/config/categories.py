# coding: utf-8

"""
Definition of categories.
"""

import order as od

from columnflow.config_util import add_category
from columnflow.util import DotDict
from columnflow.util import maybe_import
np = maybe_import("numpy")


def add_categories(config: od.Config,
                   channel = None) -> None:
    
    
    def add_base_categories(config, channel, category_map, base_selection=[]):
        base_cat = config.get_category('_'.join(('cat',channel)))      
        for i, (cat_name, cat) in enumerate(category_map.items()):
            kwargs = {
                'name'      : '_'.join((base_cat.name, cat_name)),
                'selection' : base_selection + cat.selection,
                'id'        : 100*(i+1)+ base_cat.id,
                'label'     :' '.join((base_cat.label.split(' ')[0], cat.label if 'label' in cat else cat_name))
            }
            if 'aux' in cat.keys():
                kwargs['aux'] = {}
                for (aux_spec, aux_content) in cat.aux.items():
                    if 'regs' in aux_spec:
                        kwargs['aux'][aux_spec] = {key: '_'.join((base_cat.name, val)) for (key,val) in aux_content.items()}
                    else:
                        kwargs['aux'][aux_spec] = aux_content

            add_category(config, **kwargs)
            
    
    def add_child_category(config, parent_cat, child_cat, child_name):
        max_cat_id = np.max(config.categories.ids())
        kwargs = {
                    'name'      : '__'.join((parent_cat.name, child_name)),
                    'selection' : parent_cat.selection + child_cat.selection,
                    'id'        : int(max_cat_id+1),
                    'label'     : ' '.join((parent_cat.label, child_cat.label if 'label' in child_cat else child_name))
                }
        if parent_cat.aux:
            kwargs['aux'] = {}
            for (aux_key, aux_content) in parent_cat.aux.items():
                if 'regs' in aux_key:
                    reg_map = parent_cat.aux[aux_key]
                    add_tag = lambda cat_name, tag=child_name : '__'.join((cat_name, tag))
                    reg_map_tagged = dict(zip(reg_map.keys(), map(add_tag, reg_map.values())))
                    kwargs['aux'][aux_key] = reg_map_tagged
                else:
                    kwargs['aux'][aux_key] = aux_content
        add_category(config, **kwargs)
    
    
    def create_child_categories(config, parent_categories, child_category_map):
        for cat_name in parent_categories:
            #skip 0-level categories that are used to define channelss
            if cat_name in ['incl', 'cat_mutau', 'cat_etau']: continue
            parent_cat = config.get_category(cat_name)
            for child_name, child_cat in child_category_map.items():
                add_child_category(config, parent_cat, child_cat, child_name)
    
    """
    Adds all categories to a *config*.
    ids from 1 to 9 are reserved for channels
    """
    
    add_category(
        config,
        name="incl",
        id=1,
        selection=["cat_incl"],
        label="inclusive",
    )
    if channel=='mutau':
        add_category(
            config,
            name="cat_mutau",
            id=2,
            selection=["cat_mutau"],
            label=r"$\mu\tau$ inclusive",)
        
    if channel=='etau':
        add_category(
            config,
            name="cat_etau",
            id=3,
            selection=["cat_etau"],
            label=r"$e\tau$ inclusive")

    
    
    # add_category(
    #         config,
    #         name="cat_mutau_jet_veto_check_sr",
    #         id=4,
    #         selection=["cat_mutau","jet_veto_maps_jets"],
    #         label=r"$\mu\tau$ jet veto check",)
    
    #Define initial category map with selections and call the function
    #Don't change this part: it is important for fake factor method    
    base_selection = [f'cat_{channel}','tau_eta2p3']
    
    category_map  = DotDict.wrap({
        "sr"            : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_iso", "os_charge"],
                            'label'     : "signal region",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    : 'abcd_ar',
                                               'dr_num': 'abcd_dr_num',
                                               'dr_den':  'abcd_dr_den',
                                               },
                                           #fake factor categories
                                           'ff_regs': {
                                               "ar_wj"      : "ar_wj",
                                               "dr_num_wj"  : "dr_num_wj",
                                               "dr_den_wj"  : "dr_den_wj",
                                               "ar_qcd"      : "ar_qcd",
                                               "dr_num_qcd"  : "dr_num_qcd",
                                               "dr_den_qcd"  : "dr_den_qcd",
                                               "ar_yields"   : "ar_yields",
                                               #categories for closure tests
                                               "dr_den_wj_w_ff": "dr_den_wj_w_ff",
                                               "dr_den_qcd_w_ff": "dr_den_qcd_w_ff",
                                           },},},
        #categories for jet fakes estimation via classic Fake Factor method 
        "ar_wj"          : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso", "os_charge"],
                            'aux'       : {'apply_ff': 'wj'}},
        "dr_num_wj"      : {'selection' : ['mt_inv_cut',  "deep_tau_wp",      "lep_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_wj"      : {'selection' : ['mt_inv_cut',  "deep_tau_inv_wp",  "lep_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_wj_w_ff" : {'selection' : ['mt_inv_cut',  "deep_tau_inv_wp",  "lep_iso", "tau_no_fakes", "os_charge"],
                            'aux'       : {'apply_ff': 'wj'}},
        
        "ar_qcd"         : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso"    , "os_charge"],
                            'aux'       : {'apply_ff': 'qcd'}},
        "dr_num_qcd"     : {'selection' : ['mt_cut',      "deep_tau_wp",      "lep_inv_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_qcd"     : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_inv_iso", "tau_no_fakes", "os_charge"],},
        "dr_den_qcd_w_ff": {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_inv_iso", "tau_no_fakes", "os_charge"],
                            'aux'       : {'apply_ff': 'qcd'}},
        "ar_yields"      : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso", "os_charge"],},
        #categories for QCD estimation via classic ABCD method 
        "abcd_ar"       : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "same sign region"},
        "abcd_dr_num"   : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_inv_iso", "os_charge"]},
        "abcd_dr_den"   : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_inv_iso", "ss_charge"]},
        
        "sr_no_mt"      : { 'selection' : ["deep_tau_wp", "lep_iso", "os_charge"],
                            'label'     : "signal region no mt",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    :  'abcd_ar_no_mt',
                                               'dr_num':  'abcd_dr_num_no_mt',
                                               'dr_den':  'abcd_dr_den_no_mt',
                                               },
                                           },},
        #categories for QCD estimation via classic ABCD method 
        "abcd_ar_no_mt"       : { 'selection' : ["deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "ss region no mt"},
        "abcd_dr_num_no_mt"   : { 'selection' : ["deep_tau_wp", "lep_inv_iso", "os_charge"]},
        "abcd_dr_den_no_mt"   : { 'selection' : ["deep_tau_wp", "lep_inv_iso", "ss_charge"]},
    })
    
    add_base_categories(config, channel, category_map, base_selection)
    #Add child categories to base categories
    child_category_map  = DotDict.wrap({
        
        "dm0"     : {'selection': ["pnet_dm0"], 'label': f" \n tau PNet DM = 0",},
        "dm1"     : {'selection': ["pnet_dm1"], 'label': f" \n tau PNet DM = 1",},
        #"dm2"     : {'selection': ["pnet_dm2"], 'label': f" \n tau PNet DM = 2",},
        #"dm10"    : {'selection': ["pnet_dm10"], 'label': f" \n tau PNet DM = 10",},
        #"dm11"    : {'selection': ["pnet_dm11"], 'label': f" \n tau PNet DM = 11",},
        
        "tau2pi"     : {'selection': ["pnet_dm0","tau_ip_cut"], 'label': r"$\mu \pi$",},
        "tau2rho"    : {'selection': ["pnet_dm1", "hps_dm1", "pion_E_split_cut"], 'label': r"$\mu \rho$",},

        #"nj0"    : {'selection' : ["njets_eq0"], 'label'     : f" \n $n_{{jets}}= 0$",},
        #"nj1"    : {'selection' : ["njets_eq1"], 'label'     : f" \n $n_{{jets}}= 1$",},
        #"nj2"    : {'selection' : ["njets_eq2"], 'label'     : f" \n $n_{{jets}}\geq 2$",},
        
        # "nj0_dm0"    : {'selection' : ["njets_eq0", "pnet_dm0"],
        #                      'label'     : f" \n $n_{{jets}}= 0$ \n tau PNet DM = 0",},
        
        # "nj1_dm0"    : {'selection' : ["njets_eq1", "pnet_dm0"],
        #                     'label'     : f" \n $n_{{jets}}= 1$ \n tau PNet DM = 0",},
        
        # "nj2_dm0"    : {'selection' : ["njets_geq2","pnet_dm0"],
        #                     'label'     : f" \n $n_{{jets}}\geq 2$ \n tau PNet DM = 0",},
        
        # # DM1
        # "nj0_dm1"    : {'selection' : ["njets_eq0", "pnet_dm1"],
        #                     'label'     : f" \n $n_{{jets}}= 0$ \n tau PNet DM = 1",},
        
        # "nj1_dm1"    : {'selection' : ["njets_eq1", "pnet_dm1"],
        #                     'label'     : f" \n $n_{{jets}}= 1$ \n tau PNet DM = 1",},
        
        # "nj2_dm1"    : {'selection' : ["njets_geq2","pnet_dm1"],
        #                     'label'     : f" \n $n_{{jets}}\geq 2$ \n tau PNet DM = 1",},
        
        # # DM2
        # "nj0_dm2"    : {'selection' : ["njets_eq0", "pnet_dm2"],
        #                     'label'     : f" \n $n_{{jets}}= 0$ \n tau PNet DM = 2",},
        
        # "nj1_dm2"    : {'selection' : ["njets_eq1", "pnet_dm2"],
        #                     'label'     : f" \n $n_{{jets}}= 1$ \n tau PNet DM = 2",},
        
        # "nj2_dm2"    : {'selection' : ["njets_geq2","pnet_dm2"],
        #                     'label'     : f" \n $n_{{jets}}\geq 2$ \n tau PNet DM = 2",},
        
        # # DM10
        # "nj0_dm10"    : {'selection' : ["njets_eq0", "pnet_dm10"],
        #                     'label'     : f" \n $n_{{jets}}= 0$ \n tau PNet DM = 10",},
        
        # "nj1_dm10"    : {'selection' : ["njets_eq1", "pnet_dm10"],
        #                     'label'     : f" \n $n_{{jets}}= 1$ \n tau PNet DM = 10",},
        
        # "nj2_dm10"    : {'selection' : ["njets_geq2","pnet_dm10"],
        #                     'label'     : f" \n $n_{{jets}}\geq 2$ \n tau PNet DM = 10",},
        
        # # DM11
        # "nj0_dm11"    : {'selection' : ["njets_eq0", "pnet_dm11"],
        #                     'label'     : f" \n $n_{{jets}}= 0$ \n tau PNet DM = 11",},
        
        # "nj1_dm11"    : {'selection' : ["njets_eq1", "pnet_dm11"],
        #                     'label'     : f" \n $n_{{jets}}= 1$ \n tau PNet DM = 11",},
        
        # "nj2_dm11"    : {'selection' : ["njets_geq2","pnet_dm11"],
        #                     'label'     : f" \n $n_{{jets}}\geq 2$ \n tau PNet DM = 11",},
        })
    
    create_child_categories(config,
                        parent_categories=config.categories.names(),
                        child_category_map=child_category_map)