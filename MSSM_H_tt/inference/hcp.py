# coding: utf-8

"""
Example inference model.
"""
import functools

from columnflow.inference import inference_model, ParameterType, ParameterTransformation
from columnflow.config_util import get_datasets_from_process
from MSSM_H_tt.inference.base import HCPModelBase
from MSSM_H_tt.config.mass_points import read_bdt_masses


class hcp_model(HCPModelBase):
    """
    Default statistical model for Higgs CP analysis
    """
    name = 'hcp_model'
    add_qcd = True
    processes: list = []
    config_categories: list = []
    systematics: list = []

    def init_proc_map(self) -> None:
        # mapping of process names in the datacard ("combine name") to configs and process names in a dict
        # NOTE: each value is now a list to support multiple processes per combine name
        process_vs_dataset_names = {
            "vv": [
                "www",
                "wwz",
                "zzz",
            ],
            "vvv": [
                "ww",
                "wz",
                "zz",
            ],
            "tt": [
                "tt_dl",
                "tt_fh",
                "tt_sl",
            ],
            "st": [
                "st_twchannel_tbar_fh",
                "st_twchannel_t_fh",
                "st_twchannel_tbar_dl",
                "st_twchannel_tbar_sl",
                "st_twchannel_t_dl",
                "st_twchannel_t_sl",
            ],
            "SM_higgs": [
                "h_ggf_htt",
                "h_vbf_htt",
            ],
            "vh_htt": [
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
        }

        if self.add_qcd:
            # keep QCD consistent with list semantics
            process_vs_dataset_names["QCD"] = ["qcd"]

        MASS_POINTS = read_bdt_masses()
        for m in MASS_POINTS:
            g = f"h_ggf_htt_{m}"
            b = f"bbh_htt_{m}"
            process_vs_dataset_names[g] = [g]
            process_vs_dataset_names[b] = [b]

        # store as-is (lists)
        self.proc_map = {}
        for combine_name, proc_names in process_vs_dataset_names.items():
            self.proc_map[combine_name] = proc_names

    def init_categories(self) -> None:
        ch = self.config[0].channels.names()[0]
        lep_name = ch.replace('tau', '')
        data_datasets = []
        for the_dataset in self.config[0].datasets.names():
            if f"data_{lep_name}_" in the_dataset:
                data_datasets.append(the_dataset)

        MASS_POINTS = read_bdt_masses()
        for mass in MASS_POINTS:
            self.add_category(
                f"cat_emu_sr__bdt_sig_M{mass}",
                config_category=f"cat_emu_sr__bdt_sig_M{mass}",
                config_variable=f"bdt_raw_score_sig_M{mass}",
                config_data_datasets=data_datasets,
                mc_stats=True,
                empty_bin_value=0.0,
            )

    def init_processes(self) -> None:
        config_inst = self.config[0]

        for combine_name, procs in self.proc_map.items():
            # normalize to a list
            proc_names = procs if isinstance(procs, (list, tuple, set)) else [procs]

            # QCD remains data-driven only if it's exactly ["qcd"]
            is_data_driven = (len(proc_names) == 1 and proc_names[0] == "qcd")

            is_signal = False
            dataset_names = []

            if not is_data_driven:
                for p in proc_names:
                    # resolve the process (skip and warn if not in config)
                    try:
                        proc_inst = config_inst.get_process(p)
                    except Exception:
                        print(
                            f"skipping process {p} in inference model {self.cls_name}, "
                            f"not found in config {config_inst.name}"
                        )
                        continue

                    # mark as signal if any mapped process is a signal
                    pin = proc_inst.name
                    if ("h_ggf_htt" in pin) or ("bbh_htt" in pin):
                        is_signal = True

                    # collect datasets for each mapped process
                    dsets = [
                        d.name
                        for d in get_datasets_from_process(config_inst, p, strategy="all")
                    ]
                    dataset_names.extend(dsets)

                # deduplicate while preserving order
                seen = set()
                dataset_names = [d for d in dataset_names if not (d in seen or seen.add(d))]

                if not dataset_names:
                    print(
                        f"skipping combine process {combine_name} in inference model {self.cls_name}, "
                        f"no matching datasets"
                    )
                    print(f"found in config {config_inst.name}")
                    continue

            # hand the (possibly multiple) config processes to the model layer
            self.add_process(
                name=combine_name,
                config_process=proc_names,         # list supported now
                config_mc_datasets=dataset_names,  # union of all mapped processes
                is_signal=is_signal,
                data_driven=is_data_driven,
            )

    def init_parameters(self) -> None:
        # lumi
        for config_inst in self.config:
            ckey = ''
            lumi = config_inst.x.luminosity
            for unc_name in lumi.uncertainties:
                self.add_parameter(
                    unc_name,
                    type=ParameterType.rate_gauss,
                    effect=lumi.get(names=unc_name, direction=("down", "up"), factor=True),
                    process=["*", "!QCD*"],
                    group="experiment",
                )
        # self.add_shape_parameters()


@hcp_model.inference_model
def hcp_model_no_shifts(self):
    print('Producing inference models')
    # super(hcp_model_no_shifts, self).init_func()
    hcp_model.init_func.__get__(self, self.__class__)()

    # remove all parameters that require a shift source other than nominal
    for category_name, process_name, parameter in self.iter_parameters():
        if parameter.type.is_shape or any(trafo.from_shape for trafo in parameter.transformations):
            self.remove_parameter(parameter.name, process=process_name, category=category_name)

    # repeat the cleanup
    self.init_cleanup()
