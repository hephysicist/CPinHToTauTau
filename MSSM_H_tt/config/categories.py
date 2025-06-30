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
            if cat_name in ['incl', 'cat_mutau', 'cat_etau', 'cat_tautau', 'cat_emu']: continue
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
        
        
    if channel=='tautau':
        add_category(
            config,
            name="cat_tautau",
            id=2,
            selection=["cat_tautau"],
            label=r"$\mu\tau$ inclusive",)  

    if channel=='etau':
        add_category(
            config,
            name="cat_etau",
            id=3,
            selection=["cat_etau"],
            label=r"$e\tau$ inclusive")

    if channel=='emu':
        add_category(
            config,
            name="cat_emu",
            id=4,
            selection=["cat_emu"],
            label=r"$e\mu$ inclusive")
    
    #Define initial category map with selections and call the function
    #Don't change this part: it is important for fake factor method    
    base_selection = [f'cat_{channel}']
    
    category_map  = DotDict.wrap({
        "sr"            : { 'selection' : ["lep_iso", "os_charge"],
                            'label'     : "signal region",
                            'aux'       : {
                                           #qcd estimation categories
                                           'abcd_regs' : {
                                               'ar'    : 'abcd_ar',
                                               'dr_num': 'abcd_dr_num',
                                               'dr_den':  'abcd_dr_den',
                                               },
                                           },},
        #categories for QCD estimation via classic ABCD method 
        "abcd_ar"       : { 'selection' : ["lep_iso", "ss_charge"], 'label' : "same sign region"},
        "abcd_dr_num"   : { 'selection' : ["lep_inv_iso", "os_charge"]},
        "abcd_dr_den"   : { 'selection' : ["lep_inv_iso", "ss_charge"]},

    })
    
    add_base_categories(config, channel, category_map, base_selection)
    #Add child categories to base categories
    child_category_map  = DotDict.wrap({
        # TODO: Define DZ categories here 

        "nj0"    : {'selection' : ["Zero_b_jets"], 'label'     : f" \n $n_{{jets}}= 0$",},
        "nj1"    : {'selection' : ["One_b_jets"], 'label'     : f" \n $n_{{jets}}= 1$",},
        "nj2"    : {'selection' : ["At_least_2_b_jets"], 'label'     : f" \n $n_{{jets}}\geq 2$",},       
        })
    
    create_child_categories(config,
                        parent_categories=config.categories.names(),
                        child_category_map=child_category_map)

    grand_child_category_map  = DotDict.wrap({

        "dzl"    : {'selection' : ["D_zeta_cut_low"], 'label'     : f" \n $ -35 \leq D_{{\zeta}} < -10$",},
        "dzm"    : {'selection' : ["D_zeta_cut_mid"], 'label'     : f" \n $ -10 \leq D_{{\zeta}} < 30 $",},
        "dzh"    : {'selection' : ["D_zeta_cut_high"], 'label'     : f" \n $D_{{\zeta}}\geq 30$",},       
        })
    
    create_child_categories(config,
                        parent_categories=config.categories.names(),
                        child_category_map=grand_child_category_map)
