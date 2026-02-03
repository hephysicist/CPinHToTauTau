# coding: utf-8

"""
Example inference model.
"""

import functools
import law  # for multi_match in the bin-optimization model

from columnflow.inference import inference_model, ParameterType, ParameterTransformation
from columnflow.config_util import get_datasets_from_process
from MSSM_H_tt.inference.base import HCPModelBase
from MSSM_H_tt.config.mass_points import read_bdt_masses


class hcp_model(HCPModelBase):
    """
    Default statistical model for Higgs CP analysis
    """
    name = "hcp_model"
    add_qcd = True
    processes: list = []
    config_categories: list = []
    systematics: list = []

    # -------------------------------------------------------------------------
    # process map
    # -------------------------------------------------------------------------

    def init_proc_map(self) -> None:
        # mapping of process names in the datacard ("combine name") to configs and process names in a dict
        # NOTE: each value is now a list to support multiple processes per combine name
        process_vs_dataset_names = {
            "vv": [
                "ww",
                "wz",
                "zz",
            ],
            "vvv": [
                "www",
                "wwz",
                "zzz",
            ],
            "tt": [
                "tt_dl",
                "tt_fh",
                "tt_sl",
            ],
            "st": [
                "st_tchannel_tbar",
                "st_tchannel_t",
                "st_schannel_t_lep",
                "st_schannel_tbar_lep",
                "st_twchannel_tbar_fh",
                "st_twchannel_t_fh",
                "st_twchannel_tbar_dl",
                "st_twchannel_tbar_sl",
                "st_twchannel_t_dl",
                "st_twchannel_t_sl",
            ],
            "SM_higgs": [
                "h_ggf_htt_sm_prod_sm",
                "h_vbf_htt_sm",
            ],
            "vh_htt": [
                "zh_htt_flat",
                "wph_htt_flat",
                "wmh_htt_flat",
            ],
            "wj": [
                "wj",
                "wj_1j",
                "wj_2j",
                "wj_3j",
                "wj_4j",
            ],
            "dy_tt_m50": [
                "dy_tt_m50_0j",
                "dy_tt_m50_1j",
                "dy_tt_m50_2j",
            ],
            "dy_lep": [
                "dy_lep_m10to50",
                "dy_ll_m50_0j",
                "dy_ll_m50_1j",
                "dy_ll_m50_2j",
                "dy_ll_m50",
            ],
        }

        if self.add_qcd:
            # keep QCD consistent with list semantics
            process_vs_dataset_names["QCD"] = ["qcd"]

        MASS_POINTS = read_bdt_masses()
        # temporary override, keep your original behavior if you want
        for m in MASS_POINTS:
            g = f"h_ggf_htt_{m}"
            b = f"bbh_htt_{m}"
            process_vs_dataset_names[g] = [g]
            process_vs_dataset_names[b] = [b]

        # store as-is (lists)
        self.proc_map = {}
        for combine_name, proc_names in process_vs_dataset_names.items():
            self.proc_map[combine_name] = proc_names

    # -------------------------------------------------------------------------
    # categories (new API)
    # -------------------------------------------------------------------------

    def init_categories(self) -> None:
        """
        Define categories using the new InferenceModel API:

        - use `config_data` mapping instead of `config_category`, `config_variable`,
          `config_data_datasets`.
        - build one entry per config instance via `category_config_spec`.
        """
        # get list of config instances from the base class
        config_insts = getattr(self, "config_insts", None)
        if not config_insts:
            config_insts = []
            for cfg in getattr(self, "config", []):
                if isinstance(cfg, (list, tuple, set)):
                    config_insts.extend(cfg)
                else:
                    config_insts.append(cfg)

        # use the first config to deduce the lepton flavor in the channel name
        cfg0 = config_insts[0]
        ch = cfg0.channels.names()[0]
        lep_name = ch.replace("tau", "")

        MASS_POINTS = read_bdt_masses()

        for mass in MASS_POINTS:
            config_data_ggh = {}
            config_data_bbh = {}

            for config_inst in config_insts:
                data_datasets = [
                    ds_name
                    for ds_name in config_inst.datasets.names()
                    if f"data_{lep_name}_" in ds_name
                ]

                config_data_ggh[config_inst.name] = self.category_config_spec(
                    category=f"cat_emu_sr__bdt_ggh_M{mass}",
                    variable=f"bdt_raw_score_ggh_M{mass}",
                    data_datasets=data_datasets,
                )
                config_data_bbh[config_inst.name] = self.category_config_spec(
                    category=f"cat_emu_sr__bdt_bbh_M{mass}",
                    variable=f"bdt_raw_score_bbh_M{mass}",
                    data_datasets=data_datasets,
                )

            self.add_category(
                name=f"cat_emu_sr__bdt_ggh_M{mass}",
                config_data=config_data_ggh,
                mc_stats=True,
                empty_bin_value=0.0,
            )

            self.add_category(
                name=f"cat_emu_sr__bdt_bbh_M{mass}",
                config_data=config_data_bbh,
                mc_stats=True,
                empty_bin_value=0.0,
            )

    # -------------------------------------------------------------------------
    # processes (new API – with proper QCD handling)
    # -------------------------------------------------------------------------

    def init_processes(self) -> None:
        """
        Build processes using the new `config_data` + `process_config_spec` API.

        For each combine process (e.g. "tt"):

        - for QCD (data-driven), add a process without any config_data;
        - for all others:
          * loop over all config_insts,
          * collect all MC datasets belonging to all mapped config processes,
          * create one `process_config_spec` per config_inst with `mc_datasets`
            equal to the union of all those datasets,
          * set `is_signal` if any mapped process name looks like signal.
        """
        config_insts = getattr(self, "config_insts", None)
        if not config_insts:
            config_insts = []
            for cfg in getattr(self, "config", []):
                if isinstance(cfg, (list, tuple, set)):
                    config_insts.extend(cfg)
                else:
                    config_insts.append(cfg)

        for combine_name, procs in self.proc_map.items():
            proc_names = procs if isinstance(procs, (list, tuple, set)) else [procs]

            # special case: QCD is data-driven, there is no config process "qcd"
            is_qcd = (combine_name == "QCD") or (
                len(proc_names) == 1 and proc_names[0] == "qcd"
            )
            if is_qcd:
                # create a process without any config_data; downstream tasks treat it as data-driven
                self.add_process(
                    name=combine_name,
                    is_signal=False,
                )
                continue

            # standard MC / signal processes
            is_signal = False
            config_data = {}

            for config_inst in config_insts:
                dataset_names = []

                for p in proc_names:
                    try:
                        proc_inst = config_inst.get_process(p)
                    except Exception:
                        # only warn for non-QCD processes that truly don't exist in this config
                        print(
                            f"skipping process {p} in inference model {self.cls_name}, "
                            f"not found in config {config_inst.name}"
                        )
                        continue

                    pin = proc_inst.name
                    if ("h_ggf_htt" in pin) or ("bbh_htt" in pin):
                        is_signal = True

                    dsets = [
                        d.name
                        for d in get_datasets_from_process(
                            config_inst,
                            p,
                            strategy="all",
                        )
                    ]
                    dataset_names.extend(dsets)

                # de-duplicate, keep order
                seen = set()
                dataset_names = [
                    d for d in dataset_names if not (d in seen or seen.add(d))
                ]

                if dataset_names:
                    config_data[config_inst.name] = self.process_config_spec(
                        process=None,              # rely on datasets only
                        mc_datasets=dataset_names,
                    )

            if not config_data:
                print(
                    f"skipping combine process {combine_name} in inference model {self.cls_name}, "
                    f"no matching datasets in any config"
                )
                continue

            # register the process with the new API
            self.add_process(
                name=combine_name,
                is_signal=is_signal,
                config_data=config_data,
            )

    # -------------------------------------------------------------------------
    # parameters (groups + lumi, as before)
    # -------------------------------------------------------------------------

    def init_parameters(self) -> None:
        # define common parameter groups
        if hasattr(self, "add_parameter_group"):
            for group_name in [
                "experiment",
                "theory",
                "rate_nuisances",
                "shape_nuisances",
                "signal_norm_xs",
                "signal_norm_xsbr",
            ]:
                if not self.has_parameter_group(group_name):
                    self.add_parameter_group(group_name)

        # collect config instances
        config_insts = getattr(self, "config_insts", None)
        if config_insts is None:
            config_insts = []
            for cfg in getattr(self, "config", []):
                if isinstance(cfg, (list, tuple, set)):
                    config_insts.extend(cfg)
                else:
                    config_insts.append(cfg)

        # -------------------------
        # lumi uncertainties (ADD ONCE)
        # -------------------------
        # union of all lumi nuisances across configs
        lumi_uncs = []
        seen_uncs = set()
        for cfg in config_insts:
            for unc_name in cfg.x.luminosity.uncertainties:
                if unc_name not in seen_uncs:
                    seen_uncs.add(unc_name)
                    lumi_uncs.append(unc_name)

        # process selector (for new API)
        process = ["*", "!QCD*"]
        kwargs = {}
        if hasattr(self, "process_matches"):
            # usually you want the UNION across configs; use any() as match mode
            process = self.process_matches(configs=config_insts, skip_qcd=True)
            kwargs["process_match_mode"] = any

        group = (
            ["experiment", "rate_nuisances"]
            if hasattr(self, "add_parameter_group")
            else "experiment"
        )

        # add each lumi nuisance once; sanity-check effects across configs
        for unc_name in lumi_uncs:
            ref_eff = None
            ref_cfg = None
            for cfg in config_insts:
                lumi = cfg.x.luminosity
                if unc_name not in lumi.uncertainties:
                    continue
                eff = lumi.get(names=unc_name, direction=("down", "up"), factor=True)
                if ref_eff is None:
                    ref_eff, ref_cfg = eff, cfg.name
                else:
                    if eff != ref_eff:
                        raise ValueError(
                            f"lumi nuisance '{unc_name}' has different effects across configs "
                            f"(e.g. {ref_cfg}: {ref_eff}, {cfg.name}: {eff}). "
                            "Either harmonize the lumi config or use per-config parameter names."
                        )

            # now register once
            self.add_parameter(
                unc_name,
                type=ParameterType.rate_gauss,
                effect=ref_eff,
                process=process,
                group=group,
                **kwargs,
            )


# -------------------------------------------------------------------------
# inference-model variants
# -------------------------------------------------------------------------


@hcp_model.inference_model
def hcp_model_no_shifts(self):
    """
    Analogue of `default_no_shifts` in the HH model.
    """
    print("Producing inference models without shape-based shifts")

    # initialize as in the nominal model
    super(hcp_model_no_shifts, self).init_func()

    # remove all parameters that require a shift source other than nominal
    for category_name, process_name, parameter in self.iter_parameters():
        remove = (
            (parameter.type.is_shape and not parameter.transformations.any_from_rate)
            or (parameter.type.is_rate and parameter.transformations.any_from_shape)
        )
        if remove:
            self.remove_parameter(
                parameter.name,
                process=process_name,
                category=category_name,
            )

    # repeat the cleanup
    self.init_cleanup()


@hcp_model.inference_model(empty_bin_value=0)
def hcp_model_bin_opt(self):
    """
    Bin-optimization variant, analogous to `default_bin_opt` in the HH model.
    """
    # set everything up as in the default model
    super(hcp_model_bin_opt, self).init_func()

    keep_parameters = {
        "BR_*",
        "QCDscale_*",
        "bbH_norm_*",
        "lumi_*",
        "CMS_bbtt_eff_trig_*",
        "CMS_btag_*",
        "CMS_eff_e_*",
        "CMS_eff_mu_*",
        "CMS_eff_t_*",
        "CMS_top_pT_reweighting",
        "pdf_*",
        "ps_*",
        "scale_*",
    }

    for category_name, process_name, parameter in self.iter_parameters():
        if not law.util.multi_match(parameter.name, keep_parameters):
            self.remove_parameter(
                parameter.name,
                process=process_name,
                category=category_name,
            )

    self.init_cleanup()