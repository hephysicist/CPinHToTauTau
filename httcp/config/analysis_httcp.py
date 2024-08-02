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
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    "$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

# ------------- #
# setup configs #
# ------------- #

#------------------------ Run3 2022 preEE signal samples with TauSpinner weights ----------------------- #
from httcp.config.run3_2022_preEE_tau_spinner import add_run3_2022_preEE_tau_spinner
from cmsdb.campaigns.run3_2022_preEE_tau_spinner import campaign_run3_2022_preEE_tau_spinner

add_run3_2022_preEE_tau_spinner(
    analysis_httcp,
    campaign_run3_2022_preEE_tau_spinner.copy(),
    config_name=f"{campaign_run3_2022_preEE_tau_spinner.name}_limited",
    config_id=4,
    limit_dataset_files=1)

add_run3_2022_preEE_tau_spinner(
    analysis_httcp,
    campaign_run3_2022_preEE_tau_spinner.copy(),
    config_name=f"{campaign_run3_2022_preEE_tau_spinner.name}",
    config_id=5,)
# ------------------------------------------------------------- #
