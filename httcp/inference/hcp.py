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
    
    def init_proc_map(self) -> None:
        # mapping of process names in the datacard ("combine name") to configs and process names in a dict
        
        name_map = dict([
            #("ZEE",'dy_z2ee'),
            #("ZMM",'dy_z2mumu'),
            #("ZTT",'dy_z2tautau'),
            ('ZLL','dy_lep'),
            ("W","wj"),
            ("VV",'vv'),
            ("TT","tt"),
            ('ST','st'),
            ('ggH_htt_sm', 'h_ggf_htt_sm'),
            ('ggH_htt_mm', 'h_ggf_htt_mm'),
            ('ggH_htt_cpo', 'h_ggf_htt_cpo'),
        ])
        if self.add_qcd:
            name_map["QCD"] = "qcd"
        
        self.proc_map = {}
        config_inst = self.config[0]
        for combine_name, proc_name in name_map.items():
            self.proc_map[combine_name] = proc_name

    def init_categories(self) -> None:
        ch= self.config[0].channels.names()[0]
        lep_name = ch.replace('tau','')
        data_datasets = []
        for the_dataset in self.config[0].datasets.names():
            if f"data_{lep_name}_" in the_dataset: data_datasets.append(the_dataset)
        for cat in ["tau2pi", "tau2rho"]:
            for bdt_reg in ["hig_cat_0","hig_cat_1","hig_cat_2"]:
                the_cat = self.config[0].get_category(f"cat_{ch}_sr__{cat}__{bdt_reg}")
                self.add_category(
                    f"cat_{self.campaign_key}_{ch}_sr_{cat}_{bdt_reg}",
                    config_category=f"cat_{ch}_sr__{cat}__{bdt_reg}",
                    config_variable=the_cat.x.fit_var,
                    config_data_datasets=data_datasets,
                    mc_stats = True,
                    empty_bin_value=0.0
                )
    # def init_categories(self) -> None:
    #     ch= self.config[0].channels.names()[0]
    #     lep_name = ch.replace('tau','')
    #     data_datasets = []
    #     for the_dataset in self.config[0].datasets.names():
    #         if f"data_{lep_name}_" in the_dataset: data_datasets.append(the_dataset)
    #     for cat in ["tau2pi", "tau2rho"]:
    #             the_cat = self.config[0].get_category(f"cat_{ch}_sr__{cat}__bdt_hig")
    #             self.add_category(
    #                 f"cat_{self.campaign_key}_{ch}_sr_{cat}__bdt_hig",
    #                 config_category=f"cat_{ch}_sr__{cat}__bdt_hig",
    #                 config_variable=the_cat.x.fit_var,
    #                 config_data_datasets=data_datasets,
    #                 mc_stats = True
    #             )

    def init_processes(self) -> None:
        config_inst = self.config[0]
        for combine_name, proc_name in self.proc_map.items():
            is_data_driven = (proc_name == "qcd")
            is_signal = False
            dataset_names = []
            if not is_data_driven:
                proc_inst = config_inst.get_process(proc_name)
                is_signal = ("h_ggf_htt" in proc_inst.name)
                dataset_names = [
                    dataset.name
                    for dataset in get_datasets_from_process(config_inst, proc_name, strategy="all")
                ]
                if not dataset_names:
                    print(f"skipping process {proc_name} in inference model {self.cls_name}, no matching datasets ")
                    print(f"found in config {config_inst.name}")
            self.add_process(
                name=combine_name,
                config_process=proc_name, 
                config_mc_datasets=dataset_names,
                is_signal = is_signal,
                data_driven=is_data_driven,
                
            )

    def init_parameters(self) -> None:
        # general groups
        #self.add_parameter_group("experiment")
        #self.add_parameter_group("theory")

        # lumi
        for config_inst in self.config:
            ckey = ''
            lumi = config_inst.x.luminosity
            for unc_name in lumi.uncertainties:
                self.add_parameter(
                    unc_name,
                    type=ParameterType.rate_gauss,
                    effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
                    process=[f"*", "!QCD*"],
                    group="experiment",
                )
                    

@hcp_model.inference_model
def hcp_model_no_shifts(self):
    print('Producing inference models')
    #super(hcp_model_no_shifts, self).init_func()
    hcp_model.init_func.__get__(self, self.__class__)()
    
    # remove all parameters that require a shift source other than nominal
    for category_name, process_name, parameter in self.iter_parameters():
        if parameter.type.is_shape or any(trafo.from_shape for trafo in parameter.transformations):
            self.remove_parameter(parameter.name, process=process_name, category=category_name)

    # repeat the cleanup
    self.init_cleanup()
    






























@inference_model
def example(self):

    #
    # categories
    

    self.add_category(
        "cat_mutau_sr__tau2pi",
        config_category="cat_mutau_sr__tau2pi",
        config_variable="phi_cp_mu_pi",
        config_data_datasets=["data_singlemu_C", "data_mu_C", "data_mu_D"],
        mc_stats=True,
    )

    self.add_category(
        "cat_mutau_sr",
        config_category="cat_mutau_sr",
        config_variable="phi_cp_incl",
        config_data_datasets=["data_singlemu_C", "data_mu_C", "data_mu_D"],
        mc_stats=True,
    )

    self.add_category(
        "cat_mutau_sr__tau2rho",
        config_category="cat_mutau_sr__tau2rho",
        config_variable="phi_cp_mu_rho",
        config_data_datasets=["data_singlemu_C", "data_mu_C", "data_mu_D"],
        mc_stats=True,
    )


    # TODO: think about defining a well motivated CR
    # self.add_category(
    #     "cat_mutau_abcd_ar",
    #     config_category="cat_mutau_abcd_ar",
    #     config_variable="mutau_mt",
    #     config_data_datasets=["data_singlemu_C", "data_mu_C", "data_mu_D"],
    #     mc_stats=True,
    # )

    
    # processes and datasets

    process_vs_dataset_names = {
        "data": ["data_singlemu_C", "data_mu_C", "data_mu_D"],       
    
        #Drell-Yan
        "dy_z2ee": ["dy_lep_madgraph"],
        "dy_z2mumu": ["dy_lep_madgraph"],
        "dy_z2tautau": ["dy_lep_madgraph"],
  
        "wj": ["wj_incl_madgraph"],
        #diboson
        "vv": ["ww", "wz", "zz"], #diboson inclusive
        #ttbar
        "tt": ["tt_sl","tt_dl","tt_fh"], #ttbar inclusive
        #single top
        "st": ["st_twchannel_t_sl", "st_twchannel_tbar_sl", "st_twchannel_tbar_dl", "st_tchannel_tbar", "st_tchannel_t", "st_schannel_t_lep", "st_schannel_tbar_lep"], #single top inclusive
        #signal
        #"h_ggf_htt": ["h_ggf_htt_filtered"], #SM Higgs signal SM 
        "h_ggf_htt_cpo": ["h_ggf_htt_cpo_filtered"], # CP-odd 
        "h_ggf_htt_mm": ["h_ggf_htt_mm_filtered"], # Maximal Mixing
        "h_ggf_htt_sm": ["h_ggf_htt_sm_filtered"], #CP-even 
        #data-driven backgrounds
        "qcd": [""], 
        #"jet_fakes": [""], # needs FF method
    }
 
    find_datasets = functools.partial(get_datasets_from_process, self.config_inst, strategy="all")

    for process_name, dataset_names in process_vs_dataset_names.items():
 
        is_signal = False
        data_driven = False

        if "h_ggf_htt" in process_name: 
            is_signal = True
        if process_name == "qcd" or process_name == "jet_fakes": 
            data_driven = True
            

        self.add_process(
            process_name,
            config_process=process_name, 
            config_mc_datasets=dataset_names,
            is_signal=is_signal,
            data_driven=data_driven,
        )

    #
    # parameters
    #

    # groups
 #   self.add_parameter_group("experiment")
 #   self.add_parameter_group("theory")

    # lumi
    lumi = self.config_inst.x.luminosity
    for unc_name in lumi.uncertainties:
        self.add_parameter(
            unc_name,
            type=ParameterType.rate_gauss,
            effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
            transformations=[ParameterTransformation.symmetrize],
        )


@inference_model
def example_no_shapes(self):
    # same initialization as "example" above
    example.init_func.__get__(self, self.__class__)()

    #
    # remove all shape parameters
    #

    for category_name, process_name, parameter in self.iter_parameters():
        if parameter.type.is_shape or any(trafo.from_shape for trafo in parameter.transformations):
            self.remove_parameter(parameter.name, process=process_name, category=category_name)
