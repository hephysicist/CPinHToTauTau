# coding: utf-8

"""
Configuration of the MSSM analysis.
"""
import law
import order as od

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_MSSM_H_tt_skim_2025_v1 = ana = od.Analysis(
    name="analysis_MSSM_H_tt_skim_2025_v1",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    #"$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

# ------------- #
# setup configs #
# ------------- #

from MSSM_H_tt.config.config_run3_skim_2025_v1 import add_run3
# ------------------------------------------------------------- #

channels = ['etau','mutau','emu','tautau']

#------------------------ Run3 2022 preEE samples ----------------------- #
from cmsdb.campaigns.run3_2022_preEE_nano_tau_skim_2025_v1 import campaign_run3_2022_preEE_nano_tau_skim_2025_v1
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2022_preEE_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2022_preEE_{value}_limited",
        config_id=6+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2022_preEE_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2022_preEE_{value}",
        config_id=10+counter,)
# -------------------------------------------------------------------------------------------------- #

#------------------------ Run3 2022 postEE samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2022_postEE_nano_tau_skim_2025_v1 import campaign_run3_2022_postEE_nano_tau_skim_2025_v1
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2022_postEE_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2022_postEE_{value}_limited",
        config_id=14+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2022_postEE_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2022_postEE_{value}",
        config_id=18+counter,)
# -------------------------------------------------------------------------------------------------- #


#------------------------ Run3 2023 preBPix samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2023_preBPix_nano_tau_skim_2025_v1 import campaign_run3_2023_preBPix_nano_tau_skim_2025_v1
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2023_preBPix_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2023_preBPix_{value}_limited",
        config_id=22+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2023_preBPix_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2023_preBPix_{value}",
        config_id=26+counter,)
# -------------------------------------------------------------------------------------------------- #

#------------------------ Run3 2023 postBPix samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2023_postBPix_nano_tau_skim_2025_v1 import campaign_run3_2023_postBPix_nano_tau_skim_2025_v1
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2023_postBPix_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2023_postBPix_{value}_limited",
        config_id=30+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_MSSM_H_tt_skim_2025_v1,
        campaign_run3_2023_postBPix_nano_tau_skim_2025_v1.copy(),
        channel=value,
        config_name=f"run3_2023_postBPix_{value}",
        config_id=34+counter,)
# -------------------------------------------------------------------------------------------------- #