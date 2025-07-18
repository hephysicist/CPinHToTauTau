# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import law
import order as od

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_httcp = ana = od.Analysis(
    name="analysis_httcp",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
# ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]

# default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
# if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
#     ana.x.bash_sandboxes.append(default_sandbox.name)


ana.x.bash_sandboxes = [
    "$CF_BASE/sandboxes/cf.sh",
    "$CF_BASE/sandboxes/venv_columnar.sh",
    "$HTTCP_BASE/sandboxes/venv_columnar_xgb.sh",
]

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

from httcp.config.config_run3 import add_run3
# ------------------------------------------------------------- #

channels = ['mutau','etau']

#------------------------ Run3 2022 preEE samples ----------------------- #
from cmsdb.campaigns.run3_2022_preEE_nano_tau_skim_v2 import campaign_run3_2022_preEE_nano_tau_skim_v2

for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2022_preEE_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2022_preEE_{value}_limited",
        config_id=6+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2022_preEE_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2022_preEE_{value}",
        config_id=9+counter,)
# -------------------------------------------------------------------------------------------------- #

#------------------------ Run3 2022 postEE samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2022_postEE_v2_nano_tau_v14 import campaign_run3_2022_postEE_v2_nano_tau_v14
for counter, value in enumerate(channels): 
    add_run3(
        analysis_httcp,
        campaign_run3_2022_postEE_v2_nano_tau_v14.copy(),
        channel=value,
        config_name=f"run3_2022_postEE_{value}_limited",
        config_id=12+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2022_postEE_v2_nano_tau_v14.copy(),
        channel=value,
        config_name=f"run3_2022_postEE_{value}",
        config_id=15+counter,)
# -------------------------------------------------------------------------------------------------- #


#------------------------ Run3 2023 preBPix samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2023_preBPix_nano_tau_skim_v2 import campaign_run3_2023_preBPix_nano_tau_skim_v2
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2023_preBPix_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2023_preBPix_{value}_limited",
        config_id=18+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2023_preBPix_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2023_preBPix_{value}",
        config_id=21+counter,)
# -------------------------------------------------------------------------------------------------- #

#------------------------ Run3 2023 postBPix samples ------------------------------------------------- #
from cmsdb.campaigns.run3_2023_postBPix_nano_tau_skim_v2 import campaign_run3_2023_postBPix_nano_tau_skim_v2
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2023_postBPix_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2023_postBPix_{value}_limited",
        config_id=24+counter,
        limit_dataset_files=1)
for counter, value in enumerate(channels):
    add_run3(
        analysis_httcp,
        campaign_run3_2023_postBPix_nano_tau_skim_v2.copy(),
        channel=value,
        config_name=f"run3_2023_postBPix_{value}",
        config_id=27+counter,)
# -------------------------------------------------------------------------------------------------- #

