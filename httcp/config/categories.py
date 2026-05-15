# coding: utf-8

"""
Definition of categories.
"""

import order as od
import law
from columnflow.config_util import add_category
from law.util import DotDict
from columnflow.util import maybe_import
np = maybe_import("numpy")


def add_categories(config: od.Config,
                   channel = None) -> None:
    
    
    def add_base_categories(config, channel, category_map, base_selection=[]):
        base_cats = []
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
                        if aux_spec == 'fit_var' and isinstance(aux_content, str):
                            aux_content = [aux_content]
                        kwargs['aux'][aux_spec] = aux_content
            base_cats.append(kwargs['name'])
            add_category(config, **kwargs)
        return base_cats
            
    
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
        if ('aux' in child_cat.keys()) and parent_cat.aux:
            for (aux_key, aux_content) in child_cat['aux'].items():
                if aux_key == 'fit_var' and isinstance(aux_content, str):
                    aux_content = [aux_content]
                kwargs['aux'][aux_key] = aux_content
                
        add_category(config, **kwargs)
        return kwargs['name']
    
    
    def create_child_categories(config, parent_categories, child_category_map):
        out_cats = []
        for cat_name in parent_categories:
            #skip 0-level categories that are used to define channelss
            if cat_name in ['incl', 'cat_mutau', 'cat_etau']: continue
            parent_cat = config.get_category(cat_name)
            for child_name, child_cat in child_category_map.items():
                full_name = add_child_category(config, parent_cat, child_cat, child_name)
                out_cats.append(full_name)
        return out_cats
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

    #Define initial category map with selections and call the function
    #Don't change this part: it is important for fake factor method    
    base_selection = [f'cat_{channel}','tau_eta2p3','muon_ip_cut']
    
    category_map  = DotDict.wrap({
        "sr"            : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_iso", "os_charge"],
                            'label'     : r"incl. SR",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    : 'abcd_ar',
                                               #'dr_num': 'abcd_dr_num',
                                               #'dr_den':  'abcd_dr_den',
                                               },
                                           #fake factor categories
                                           'ff_regs': {
                                               "ar_wj"       : "ar_wj",
                                               "dr_num_wj"   : "dr_num_wj",
                                               "dr_den_wj"   : "dr_den_wj",
                                               "ar_qcd"      : "ar_qcd",
                                               "dr_num_qcd"  : "dr_num_qcd",
                                               "dr_den_qcd"  : "dr_den_qcd",
                                               "ar_yields"   : "ar_yields",
                                                #categories for closure tests
                                                "dr_den_wj_w_ff": "dr_den_wj_w_ff",
                                                #"dr_den_wj_w_ff_fakes":"dr_den_wj_w_ff_fakes",
                                                "dr_den_qcd_w_ff": "dr_den_qcd_w_ff",
                                                #"dr_den_qcd_w_ff_fakes":"dr_den_qcd_w_ff_fakes",
                                                
                                           },},},
        #categories for jet fakes estimation via classic Fake Factor method 
       
        "dr_num_wj"      : {'selection' : ['mt_inv_cut',  "deep_tau_wp",      "lep_iso", "ss_charge"],},
        "dr_den_wj"      : {'selection' : ['mt_inv_cut',  "deep_tau_inv_wp",  "lep_iso", "ss_charge"],},
        "dr_den_wj"      : {'selection' : ['mt_inv_cut',  "deep_tau_inv_wp",  "lep_iso", "ss_charge"],},
        "dr_den_wj_w_ff"       : {'selection' : ['mt_inv_cut',  "deep_tau_inv_wp",  "lep_iso", "ss_charge"], 'aux':{'apply_ff': 'wj'}},
        
      
        "dr_num_qcd"            : {'selection' : ['mt_cut',      "deep_tau_wp",      "lep_inv_iso", "ss_charge"],},
        "dr_den_qcd"            : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_inv_iso", "ss_charge"],},
        "dr_den_qcd_w_ff"       : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_inv_iso", "ss_charge"], 'aux': {'apply_ff': 'qcd'}},
       
        "ar_yields"             : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso",     "os_charge"],},
        "ar_wj"          : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso", "os_charge"], 'aux': {'apply_ff': 'wj'}},
        "ar_qcd"                : {'selection' : ['mt_cut',      "deep_tau_inv_wp",  "lep_iso",     "os_charge"], 'aux': {'apply_ff': 'qcd'}},
        
        # #categories for QCD estimation via classic ABCD method 
        "abcd_ar"       : { 'selection' : ['mt_cut', "deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "same sign region"},
        
        "sr_no_mt"      : { 'selection' : ["deep_tau_wp", "lep_iso", "os_charge"],
                            'label'     : "inclusive, no mT",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    :  'abcd_ar_no_mt',},
                                           },},
        #categories for QCD estimation via classic ABCD method 
        "abcd_ar_no_mt"       : { 'selection' : ["deep_tau_wp", "lep_iso", "ss_charge"], 'label' : "same sign, no mT"},
    })
    
    base_cats = add_base_categories(config, channel, category_map, base_selection)
    #Add child categories to base categories
    
    jet_fake_map  = DotDict.wrap({
        "jet_fakes" : {'selection'  : ["tau_jet_fakes"],    'label': f"\njet fakes",},
        "prompt"    : {'selection'   : ["tau_no_jet_fakes"], 'label': f"\nprompt lep."},
        "fake_incl" : {'selection'   : [], 'label': f"\nfake incl."},
        })
    
    bdt_map  = DotDict.wrap({
        "hig"      : {'selection'  : ["bdt_cat_higgs"], 'label': f"\nbdt cat Higgs",},
        "gtau"     : {'selection'  : ["bdt_cat_gtau"], 
                      'label'      : f"\nbdt cat genuine tau",
                      'aux'        : {'fit_var': 'bdt_raw_score_4bins_gtau'},},
        "fake"     : {'selection'  : ["bdt_cat_fake"],
                      'label'      : f"\nbdt cat fakes",
                      'aux'        : {'fit_var': 'bdt_raw_score_4bins_fake'},},
        "bdt_incl" : {'selection'  : [],
                      'label': f"\nBDT incl."}, 
        })

    hig_bins_map  = DotDict.wrap({
        "cat0"    : {'selection': ["hig_cat_0"], 'label': f"\nD_H in (0.33,0.5]",},
        "cat1"    : {'selection': ["hig_cat_1"], 'label': f"\nD_H in (0.5,0.7]",},
        "cat2"    : {'selection': ["hig_cat_2"], 'label': f"\nD_H in (0.7,1.0]",},
        "cat_incl": {'selection': [], 'label': f"",},
        })
    
    tau_decays_map  = DotDict.wrap({
        "tau2pi"     : {'selection': ["pnet_dm0","tau_ip_cut"], 
                        'label': f"\n mu pi",
                        'aux'       : {'fit_var': 'phi_cp_mu_pi'},
                        },

        "tau2rho"    : {'selection': ["pnet_dm1", "hps_dm1", "pion_E_split_cut", "tau_has_em_prods"],
                        'label': f"\n mu rho",
                        'aux'       : {'fit_var': 'phi_cp_mu_rho'},
                        },

        "tau2a1"    : {'selection': ["pnet_dm2", "hps_dm1","pion_E_split_cut","tau_has_em_prods"],
                        'label': f"\n mu a1 1pr",
                        'aux'       : {'fit_var': 'phi_cp_mu_a1_1pr'},
                        },

        "tau2a1_3pr": {'selection': ["pnet_dm10","has_refit_sv"],
                        'label': f"\n mu a1 3pr",
                        'aux'       : {'fit_var': 'phi_cp_mu_a1_3pr_pv_gef'},
                        },
        })
    
    splitted_by_jet_fake = create_child_categories(config,
                                                 parent_categories=base_cats,
                                                 child_category_map=jet_fake_map)
    splitted_by_bdt = create_child_categories(config,
                                              parent_categories=splitted_by_jet_fake,
                                              child_category_map=bdt_map)
    
    sr_cats = [c for c in splitted_by_bdt if 'hig' in c]
    
    splitted_by_hig_bins = create_child_categories(config,
                                                   parent_categories=sr_cats,
                                                   child_category_map=hig_bins_map)
    
    splitted_by_dm_hig = create_child_categories(config,
                                                 parent_categories=splitted_by_hig_bins,
                                                 child_category_map=tau_decays_map)
    bdt_incl_cats = [c for c in splitted_by_bdt if 'bdt_incl' in c]
    
    splitted_by_dm_incl = create_child_categories(config,
                                                 parent_categories=bdt_incl_cats,
                                                 child_category_map=tau_decays_map)
    
    # no_split_by_dm = [c for c in config.categories.names() if ('tau2' not in c) and (c != 'cat_mutau_sr') and (c != 'incl')]
    
    # for the_name in no_split_by_dm:
    #     the_cat = config.get_category(the_name)
    #     config.remove_category(the_cat)
    