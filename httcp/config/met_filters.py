# coding: utf-8

"""
Definition of MET filter flags.
"""

import order as od

from law.util import DotDict


def add_met_filters(config: od.Config) -> None:
    """
    Adds all MET filters to a *config*.

    Resources:
    https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=157#UL_data
    """
    #https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    filters = [
        "Flag.goodVertices", # +
        "Flag.globalSuperTightHalo2016Filter", # +
        "Flag.EcalDeadCellTriggerPrimitiveFilter", # +
        "Flag.BadPFMuonFilter", # +
        "Flag.BadPFMuonDzFilter", # +
        "Flag.hfNoisyHitsFilter", # + 
        "Flag.eeBadScFilter", # +
        #"Flag.ecalBadCalibFilter", # + (not in the code of IC)
    ]
    # same filter for mc and data, but still separate
    filters = {
        "mc": filters,
        "data": filters,
    }
    config.x.met_filters = DotDict.wrap(filters)
