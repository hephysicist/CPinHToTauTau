"""
Inference base models with common functionality.
"""

import re
import abc

import order as od

from columnflow.inference import InferenceModel


class HCPModelBase(InferenceModel):
    """
    Base class for statistical models of the Higgs CP analysis inherited from Hh2bbtautau repo:
    https://github.com/uhh-cms/hh2bbtautau/blob/master/hbt/inference/base.py
    """

    def __init__(self, *args, **kwargs) -> None:
        # members that are set in init_objects
     
        self.config = args
        self.single_config: bool
        self.campaign_keys: dict[od.Config, str] = {}
        self.campaign_key: str
        self.proc_map: dict[str, dict[od.Config, str]] = {}
        super().__init__(*args, **kwargs)
        
        
    @abc.abstractmethod
    def init_proc_map(self) -> None:
        # should setup self.proc_map
        ...

    @abc.abstractmethod
    def init_categories(self) -> None:
        # should setup inference model cateogries
        ...

    @abc.abstractmethod
    def init_processes(self) -> None:
        # should setup inference model processes
        ...

    @abc.abstractmethod
    def init_parameters(self) -> None:
        # should setup inference model parameters
        ...

    def init_func(self) -> None:
        # the default initialization is split into logical parts
        self.init_objects()
        self.init_categories()
        self.init_processes()
        self.init_parameters()
        self.init_cleanup()

    def inject_era(self, config_inst: od.Config, combine_name: str) -> str:
        # helper to inject era info into combine process names
        tag = config_inst.x.tag
        year = config_inst.x.year
        return f"{combine_name}_{year}{tag}"

    def init_objects(self) -> None:
        # gather campaign identifier keys per config
        #self.single_config = len(self.config) == 1
        #for config_inst in self.config:
        #    year2 = config_inst.campaign.x.year % 100
        #    self.campaign_keys[config_inst] = f"{year2}{config_inst.campaign.x.postfix}"

        # overall campaign key
        self.campaign_key = "_".join(self.campaign_keys.values())

        # setup the process_map
        self.init_proc_map()

    def init_cleanup(self) -> None:
        self.cleanup(keep_parameters="")
        
    
    def format_systematics(self):
        """ Function to format the systematics and prepare the processes_per_* dictionaries
        """
        years = {config_inst.campaign.x.year for config_inst in self.config_insts}

        systematics_formatted = sorted({
            syst.format(year=year, bjet_cat=bjet_cat)
            for syst in self.systematics
            for year in years
            for bjet_cat in self.bjet_cats
        })

        available_procs = set(self.processes)

        self.processes_per_shape = {
            unc_formatted: available_procs
            for year in years
            for bjet_cat in self.bjet_cats
            for unc, procs in const.processes_per_shape.items()
            if (
                (unc_formatted := unc.format(year=year, bjet_cat=bjet_cat))
                in systematics_formatted and
                any(available_procs := [proc for proc in procs if proc in ["all", *self.processes]])
            )
        }
        # check that all systematics are considered and warn if not
        all_systs = set()
        for processes_per_syst in (
            self.processes_per_QCDScale,
            self.processes_per_pdf_rate,
            self.processes_per_rate_unconstrained,
            self.processes_per_shape,
        ):
            for syst, procs in processes_per_syst.items():
                if procs:
                    all_systs.add(syst)
                else:
                    logger.warning(
                        f"No processes found for systematic {syst}. "
                        "Please check your configuration.",
                    )
                    processes_per_syst.pop(syst)
        if not_considered := set(systematics_formatted) - all_systs:
            logger.warning(
                f"The following systematics were not considered in the inference model: "
                f"{', '.join(not_considered)}. Please check your configuration.",
            )
