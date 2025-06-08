# coding: utf-8

"""
Example inference model.
"""
import functools 

from columnflow.inference import inference_model, ParameterType, ParameterTransformation
from columnflow.config_util import get_datasets_from_process 



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
