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
    all_datasets = self.config_inst.campaign.datasets.names()
    data_emu    = []
    mc_datasets = []

    for sample in all_datasets:
        if ("data_egamma" in sample) or ("data_mu_" in sample) or ("data_singlemu" in sample):
            data_emu.append(sample)
        else:
            mc_datasets.append(sample)

    for name in self.config_inst.categories.names():
        if name == "cat_emu_sr__bdt_sig":
            self.add_category(
                name,
                config_category=name,
                config_variable="bdt_raw_score_sig",
                config_data_datasets=data_emu,
                mc_stats=True,
            )
                
    # processes and datasets

    process_vs_dataset_names = {
        "data": data_emu,       
        "vv": [
            "ww",
            "wz",
            "zz",
            ],
        "vvv": [
            "www",
            "wwz",
            "zzz"
            ],
        "tt": [
            "tt_dl",
            "tt_fh",
            "tt_sl",
            ],
        "st": [
            "TbarWplusto4Q",
            "TWminusto4Q",
            "TbarWplusto2L2Nu", 
            "TbarWplustoLNu2Q",
            "TWminusto2L2Nu", 
            "TWminustoLNu2Q", 
            ],
        "SM_higgs": [
            "h_ggf_htt",
            "h_vbf_htt",
            ],
        "vh_htt": [
            "vh_htt",
            "zh_htt",
            "wh_htt",
            ],
        "w_lnu": [         
            "w_lnu",     
            ],
        "dy_tautau": [
            "dy_tautau_m50toinf_0j",
            "dy_tautau_m50toinf_1j",
            "dy_tautau_m50toinf_2j",
            ],
        "dy_ll": [
            "dy_m10to50",
            "dy_m50toinf_0j",
            "dy_m50toinf_1j",
            "dy_m50toinf_2j",
            "dy_m50toinf",  
            ],
       "h_ggf_htt" : ["h_ggf_htt"],
       "bbh_htt" : ["bbh_htt"],
    }
    
    signal_masses = [60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1400, 1600, 1800, 2000, 2300, 2600, 2900, 3200, 3500]
    for m in signal_masses:
        key = f"h_ggf_htt_{m}"
        process_vs_dataset_names[key] = [key]
        key = f"bbh_htt_{m}"
        process_vs_dataset_names[key] = [key]
        
    process_vs_dataset_names["qcd"] = [""]

    #find_datasets = functools.partial(get_datasets_from_process, self.config_inst, strategy="all")

    for process_name, dataset_names in process_vs_dataset_names.items():

        is_signal = False
        data_driven = False

        if "h_ggf_htt_" in process_name: 
            is_signal = True
        if "bbh_htt_" in process_name: 
            is_signal = True
        if process_name == "qcd": #or process_name == "jet_fakes": 
            data_driven = True
            
        self.add_process(
            process_name,
            config_process=process_name, 
            config_mc_datasets=dataset_names,
            is_signal=is_signal,
            data_driven=data_driven,
        )

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
        
