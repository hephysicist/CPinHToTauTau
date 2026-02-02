# coding: utf-8

"""
Example inference model.
"""
import functools 

from columnflow.inference import inference_model, ParameterType, ParameterTransformation
from columnflow.config_util import get_datasets_from_process 
from httcp.inference.base import HCPModelBase


class hcp_model(HCPModelBase):
    """
    Default statistical model for Higgs CP analysis
    """
    
    name = 'hcp_model'
    add_qcd = True
    add_fakes = False
    processes: list = []
    config_categories: list = []
    systematics: list = []
    variable = "phi_cp_mu_a1_3pr_pv_gef"
      
    def init_proc_map(self) -> None:
        # mapping of process names in the datacard ("combine name") to configs and process names in a dict
        
        name_map = dict([
            ("ZL",'dy_ll_m50'),
            ("ZTT",'dy_tt_m50'),
            #Diboson + single top
            ("VVT",'vvt'), #use hist_hooks to calculate his process
            #ttbar
            ("TTT","tt"),
            #wjets
            ("WJ",'w'),
            #VBF Signal
            ('qqH_sm_htt125',   'h_vbf_htt_sm'),
            ('qqH_mm_htt125',   'h_vbf_htt_mm'),
            ('qqH_ps_htt125',   'h_vbf_htt_cpo'),
            #('qqH_flat_htt125', 'h_vbf_htt_flat'),
            #ggF Signal
            ('ggH_sm_prod_sm_htt125', 'h_ggf_htt_sm_prod_sm'),
            ('ggH_mm_prod_sm_htt125', 'h_ggf_htt_mm_prod_sm'),
            ('ggH_ps_prod_sm_htt125', 'h_ggf_htt_cpo_prod_sm'),
            #('ggH_flat_prod_sm_htt125', 'h_ggf_htt_flat_prod_sm'),

            #('ggH_sm_prod_mm_htt125', 'h_ggf_htt_sm_prod_mm'),
            #('ggH_mm_prod_mm_htt125', 'h_ggf_htt_mm_prod_mm'),
            #('ggH_ps_prod_mm_htt125', 'h_ggf_htt_cpo_prod_mm'),
            #('ggH_flat_prod_mm_htt125', 'h_ggf_htt_flat_prod_mm'),

            #('ggH_sm_prod_cpo_htt125', 'h_ggf_htt_sm_prod_cpo'),
            #('ggH_mm_prod_cpo_htt125', 'h_ggf_htt_mm_prod_cpo'),
            #('ggH_ps_prod_cpo_htt125', 'h_ggf_htt_cpo_prod_cpo'),
            #('ggH_flat_prod_cpo_htt125', 'h_ggf_htt_flat_prod_cpo'),

            #ZH Signal
            #('ZH_sm_htt125',   'zh_htt_sm'),
            #('ZH_mm_htt125',   'zh_htt_mm'),
            #('ZH_ps_htt125',   'zh_htt_cpo'),
            #('ZH_flat_htt125', 'zh_htt_flat'),

            #WH Signal
            #('WH_sm_htt125',   'wh_htt_sm'),
            #('WH_mm_htt125',   'wh_htt_mm'),
            #('WH_ps_htt125',   'wh_htt_cpo'),
            #('WH_flat_htt125', 'wh_htt_flat'),
        ])
        if self.add_qcd:
            name_map["QCD"] = "qcd"
        if self.add_fakes:
            name_map["JetFakes"] = "qcd"
        
        self.proc_map = {}

        for combine_name, proc_name in name_map.items():
            self.proc_map[combine_name] = proc_name

    def init_categories(self) -> None:
        cfg = self.config_insts[0]
        ch = cfg.channels.names()[0]
        lep_name = ch.replace('tau','')
        data_datasets = []
        for the_dataset in cfg.datasets.names():
            if ("data_mu_" in the_dataset) or ("data_singlemu_") in the_dataset: 
                data_datasets.append(the_dataset)

        ch_names = {
            "tau2pi":  'mupi',
            "tau2rho": 'murho',
            "tau2a1":  'mua11pr',
            "tau2a1_3pr":'mua1'
        }
        #for cat in ["tau2pi", "tau2rho", "tau2a1", "tau2a1_3pr"]:
        for cat in ["tau2a1_3pr"]:
            for bdt_reg in ["cat0","cat1","cat2"]:
                the_cat = cfg.get_category(f"cat_{ch}_sr__fake_incl__hig__{bdt_reg}__{cat}")
                ch_name = ch_names[cat]
                self.add_category(
                    f"mt_mva_higgs_{bdt_reg}_{ch_name}",
                    config_data={
                        cfg.name: self.category_config_spec(
                            category=f"cat_{ch}_sr__fake_incl__hig__{bdt_reg}__{cat}",
                            variable=the_cat.x.fit_var[0],
                            data_datasets=data_datasets)
                    },
                    mc_stats = True,
                    empty_bin_value=0.0, 
                )
        #Adding background categories
        for cat_name in ['gtau','fake']:
            the_cat = cfg.get_category(f"cat_mutau_sr__fake_incl__{cat_name}")
            if cat_name == 'gtau': dc_name = 'tau'
            else : dc_name = cat_name
            self.add_category(
                    f"mt_mva_{dc_name}",
                    config_data={
                        cfg.name: self.category_config_spec(
                            category=f"cat_mutau_sr__fake_incl__{cat_name}",
                            variable=the_cat.x.fit_var[0],
                            data_datasets=data_datasets,)
                    },
                    mc_stats = True,
                    empty_bin_value=0.0
                )

    def init_processes(self) -> None:
        cfg = self.config_insts[0]
        for combine_name, proc_name in self.proc_map.items():
            is_data_driven = (proc_name == "qcd")
            is_signal = False
            dataset_names = []
            if not is_data_driven:
                proc_inst = cfg.get_process(proc_name)
                is_signal = (("h_ggf_htt" in proc_inst.name) or 
                             ("h_vbf_htt" in proc_inst.name) or 
                             ("zh_htt" in proc_inst.name) or 
                             ("wh_htt" in proc_inst.name))
                dataset_names = [
                    dataset.name
                    for dataset in get_datasets_from_process(cfg, proc_name, strategy="all")
                ]
                if not dataset_names:
                    print(f"skipping process {proc_name} in inference model {self.cls_name}, no matching datasets ")
                    print(f"found in config {cfg.name}")
            self.add_process(
                name=combine_name,
                config_data={
                        cfg.name: self.process_config_spec(
                            process=proc_name,
                            mc_datasets=dataset_names,
                            
                        )},
                is_dynamic = is_data_driven,
                is_signal = is_signal,
            )
    def init_parameters(self) -> None:
        cfg = self.config_insts[0]
        # general groups
        self.add_parameter_group("experiment")
        self.add_parameter_group("theory")

        # lumi
        ckey = ''
        lumi = cfg.x.luminosity
        for unc_name in lumi.uncertainties:
            self.add_parameter(
                unc_name,
                type=ParameterType.rate_gauss,
                effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
                process=[f"*", "!QCD*"],
                group="experiment",
            )
        #from IPython import embed; embed()
            #self.add_shape_parameters()
                
    # def add_shape_parameters(self: InferenceModel):
    #     """
    #     Function that adds all rate parameters to the inference model
    #     """
    #     processes_per_shape = {
    #         "tau": ["all"],
    #     }
        
    #     source_per_shape = {shape: shape for shape in processes_per_shape.keys()}
        
    #     for shape_uncertainty, shape_processes in self.processes_per_shape.items():

    #         # If "all" is included, takes all processes except for the ones specified (starting with !)
    #         if "all" in shape_processes:
    #             _remove_processes = {proc[:1] for proc in shape_processes if proc.startswith("!")}
    #             shape_processes = set(self.processes) - _remove_processes

    #         param_kwargs = {
    #             "type": ParameterType.shape,
    #             "process": [self.inf_proc(proc) for proc in shape_processes],
    #         }
           
    #         shift_source = const.source_per_shape.get(shape_uncertainty, shape_uncertainty)
            
    #         param_kwargs["config_data"] = {
    #             config_inst.name: self.parameter_config_spec(
    #                 shift_source=shift_source,
    #             )
    #             for config_inst in self.config_insts
    #             if config_inst.has_shift(f"{shift_source}_up")
    #         }

    #         self.add_parameter(
    #             shape_uncertainty,
    #             **param_kwargs,
    #         )

            
    #         if "pdf" in shape_uncertainty:
    #             self.add_parameter_to_group(shape_uncertainty, "theory")
    #         else:
    #             self.add_parameter_to_group(shape_uncertainty, "experiment")
                
   

@hcp_model.inference_model
def hcp_model(self):
    print('Producing inference models')
    #super(hcp_model, self).init_func()
    hcp_model.init_func.__get__(self, self.__class__)()
    # remove all parameters that require a shift source other than nominal
    for category_name, process_name, parameter in self.iter_parameters():
        remove = (
            (parameter.type.is_shape and not parameter.transformations.any_from_rate) or
            (parameter.type.is_rate and parameter.transformations.any_from_shape)
        )
        if remove:
            self.remove_parameter(parameter.name, process=process_name, category=category_name)


    # repeat the cleanup
    self.init_cleanup()
 