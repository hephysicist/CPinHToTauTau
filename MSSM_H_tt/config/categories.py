# coding: utf-8

"""
Definition of categories.
"""

import order as od

from columnflow.config_util import add_category


def add_categories(config: od.Config,
                   channel = None) -> None:
    """
    Adds all categories to a *config*.
    ids from 1 to 9 are reserved for channels
    """
    
    
    add_category(
        config,
        name="incl",
        id=1,
        selection="cat_incl",
        label="inclusive",
    )
     
    ##############################################
    ### Main categories for the three channels ###
    ##############################################
    if channel=='mutau':
        mutau = add_category(
            config,
            name="cat_mutau",
            id=2,
            selection="cat_mutau",
            label="$\\mu\\tau$ incl",
        )
        #########################################
        ### mu-tau channel categories 0 b jets ###
        #########################################
        
        #################################
        ### SIGNAL REGION with mT cut ###
        #################################
        
        mutau_signal_reg_0_bjets = add_category(
            config,
            name="mutau_signal_reg_0_bjets",
            id=130 + mutau.id,
            selection=["cat_mutau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets",
                       ],
            label="$\\mu\\tau$ SR\n0 bjets",
            aux={'control_reg': "mutau_control_reg_0_bjets"}
        )
        mutau_control_reg_0_bjets = add_category(
            config,
            name="mutau_control_reg_0_bjets",
            id=131 + mutau.id,
            selection=["cat_mutau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       "Zero_b_jets" ,
                       ],
            label="$\\mu\\tau$ CR\n0 bjets",
        )
        mutau_signal_reg_0_bjets_endcap_tau = add_category(
            config,
            name="mutau_signal_reg_0_bjets_endcap_tau",
            id=132 + mutau.id,
            selection=["cat_mutau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
            label = "$\\mu\\tau$ SR\n0 bjets\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_0_bjets_endcap_tau"}
        )
        mutau_control_reg_0_bjets_endcap_tau = add_category(
            config,
            name="mutau_control_reg_0_bjets_endcap_tau",
            id=133 + mutau.id,
            selection=["cat_mutau"   ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
           label = "$\\mu\\tau$ SR\n0 bjets\n$\\eta_{\\tau} > 1.2$",
        )
        
        mutau_signal_reg_0_bjets_barrel_tau = add_category(
            config,
            name="mutau_signal_reg_0_bjets_barrel_tau",
            id=134 + mutau.id,
            selection=["cat_mutau"    ,
                        "os_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ SR\n0 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_0_bjets_barrel_tau"}
        )
        
        mutau_control_reg_0_bjets_barrel_tau = add_category(
            config,
            name="mutau_control_reg_0_bjets_barrel_tau",
            id=135 + mutau.id,
            selection=["cat_mutau"  ,
                        "ss_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ CR\n0 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )
        
        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################
        
        mutau_signal_reg_0_bjets_no_mt = add_category(
            config,
            name="mutau_signal_reg_0_bjets_no_mt",
            id=136 + mutau.id,
            selection=["cat_mutau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                      "Zero_b_jets" ,
                       ],
            label="$\\mu\\tau$\n0 bjets\nno $m_{T}$",
            aux={'control_reg': "mutau_control_reg_0_bjets_no_mt"}
        )
        
        mutau_control_reg_0_bjets_no_mt = add_category(
            config,
            name="mutau_control_reg_0_bjets_no_mt",
            id=137 + mutau.id,
            selection=["cat_mutau"    ,
                       "ss_charge"   ,
                       "deep_tau_wp" ,
                       "Zero_b_jets" ,
                       ],
            label="$\\mu\\tau$ CR\n0 bjets\nno $m_{T}$",
        )
        mutau_signal_reg_0_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_0_bjets_endcap_tau_no_mt",
            id=138 + mutau.id,
            selection=["cat_mutau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
            label="$\\mu\\tau$ SR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_0_bjets_endcap_tau_no_mt"}
        )
        mutau_control_reg_0_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_0_bjets_endcap_tau_no_mt",
            id=139 + mutau.id,
            selection=["cat_mutau"   ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
            label="$\\mu\\tau$ CR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )
        mutau_signal_reg_0_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_0_bjets_barrel_tau_no_mt",
            id=140 + mutau.id,
            selection=["cat_mutau"    ,
                        "os_charge"  ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ SR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_0_bjets_barrel_tau_no_mt"}
        )
        
        mutau_control_reg_0_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_0_bjets_barrel_tau_no_mt",
            id=141 + mutau.id,
            selection=["cat_mutau"  ,
                        "ss_charge"  ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ CR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
        
        #########################################
        ### mu-tau channel categories 1 b jets ###
        #########################################

        #################################
        ### SIGNAL REGION with mT cut ###
        #################################

        mutau_signal_reg_1_bjets = add_category(
            config,
            name="mutau_signal_reg_1_bjets",
            id=142 + mutau.id,
            selection=["cat_mutau"   ,
                    "os_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets",
                    ],
            label="$\\mu\\tau$ SR\n1 bjet",
            aux={'control_reg': "mutau_control_reg_1_bjets"}
        )

        mutau_control_reg_1_bjets = add_category(
            config,
            name="mutau_control_reg_1_bjets",
            id=143 + mutau.id,
            selection=["cat_mutau"    ,
                    "ss_charge"   ,
                    "mt_cut"      ,
                    "deep_tau_wp" ,
                    "One_b_jets"  ,
                    ],
            label="$\\mu\\tau$ CR\n1 bjet",
        )

        mutau_signal_reg_1_bjets_endcap_tau = add_category(
            config,
            name="mutau_signal_reg_1_bjets_endcap_tau",
            id=144 + mutau.id,
            selection=["cat_mutau"   ,
                    "os_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$\\mu\\tau$ SR\n1 bjets\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_1_bjets_endcap_tau"}
        )

        mutau_control_reg_1_bjets_endcap_tau = add_category(
            config,
            name="mutau_control_reg_1_bjets_endcap_tau",
            id=145 + mutau.id,
            selection=["cat_mutau"   ,
                    "ss_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$\\mu\\tau$ CR\n1 bjets\n$\\eta_{\\tau} > 1.2$",
        )

        mutau_signal_reg_1_bjets_barrel_tau = add_category(
            config,
            name="mutau_signal_reg_1_bjets_barrel_tau",
            id=146 + mutau.id,
            selection=["cat_mutau"    ,
                        "os_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ SR\n1 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_1_bjets_barrel_tau"}
        )

        mutau_control_reg_1_bjets_barrel_tau = add_category(
            config,
            name="mutau_control_reg_1_bjets_barrel_tau",
            id=147 + mutau.id,
            selection=["cat_mutau"  ,
                        "ss_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ CR\n1 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )

        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################

        mutau_signal_reg_1_bjets_no_mt = add_category(
            config,
            name="mutau_signal_reg_1_bjets_no_mt",
            id=148 + mutau.id,
            selection=["cat_mutau"   ,
                    "os_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    ],
            label="$\\mu\\tau$ SR\n1 bjets\nno $m_{T}$",
            aux={'control_reg': "mutau_control_reg_1_bjets_no_mt"}
        )

        mutau_control_reg_1_bjets_no_mt = add_category(
            config,
            name="mutau_control_reg_1_bjets_no_mt",
            id=149 + mutau.id,
            selection=["cat_mutau"    ,
                    "ss_charge"   ,
                    "deep_tau_wp" ,
                    "One_b_jets"  ,
                    ],
            label="$\\mu\\tau$ CR\n1 bjets\nno $m_{T}$",
        )

        mutau_signal_reg_1_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_1_bjets_endcap_tau_no_mt",
            id=150 + mutau.id,
            selection=["cat_mutau"   ,
                    "os_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$\\mu\\tau$ SR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_1_bjets_endcap_tau_no_mt"}
        )

        mutau_control_reg_1_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_1_bjets_endcap_tau_no_mt",
            id=151 + mutau.id,
            selection=["cat_mutau"   ,
                    "ss_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$\\mu\\tau$ CR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )

        mutau_signal_reg_1_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_1_bjets_barrel_tau_no_mt",
            id=152 + mutau.id,
            selection=["cat_mutau"    ,
                        "os_charge"  ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ SR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_1_bjets_barrel_tau_no_mt"}
        )

        mutau_control_reg_1_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_1_bjets_barrel_tau_no_mt",
            id=153 + mutau.id,
            selection=["cat_mutau"  ,
                        "ss_charge"  ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$\\mu\\tau$ CR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
            
        #########################################
        ### mu-tau channel categories >= 2 bjets ###
        #########################################

        #################################
        ### SIGNAL REGION with mT cut ###
        #################################

        mutau_signal_reg_2_bjets = add_category(
            config,
            name="mutau_signal_reg_2_bjets",
            id=154 + mutau.id,
            selection=["cat_mutau"        ,
                    "os_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$\\mu\\tau$ SR\n$\\\geq$ 2 bjets",
            aux={'control_reg': "mutau_control_reg_2_bjets"}
        )

        mutau_control_reg_2_bjets = add_category(
            config,
            name="mutau_control_reg_2_bjets",
            id=155 + mutau.id,
            selection=["cat_mutau"        ,
                    "ss_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$\\mu\\tau$ CR\n$\\\geq$ 2 bjets",
        )

        mutau_signal_reg_2_bjets_endcap_tau = add_category(
            config,
            name="mutau_signal_reg_2_bjets_endcap_tau",
            id=156 + mutau.id,
            selection=["cat_mutau"        ,
                    "os_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$\\mu\\tau$ SR\n$\\\geq$ 2 bjets\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_2_bjets_endcap_tau"}
        )

        mutau_control_reg_2_bjets_endcap_tau = add_category(
            config,
            name="mutau_control_reg_2_bjets_endcap_tau",
            id=157 + mutau.id,
            selection=["cat_mutau"        ,
                    "ss_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$\\mu\\tau$ CR\n$\\\geq$ 2 bjets\n$\\eta_{\\tau} > 1.2$",
        )

        mutau_signal_reg_2_bjets_barrel_tau = add_category(
            config,
            name="mutau_signal_reg_2_bjets_barrel_tau",
            id=158 + mutau.id,
            selection=["cat_mutau"        ,
                        "os_charge"      ,
                        "mt_cut"         ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$\\mu\\tau$ SR\n$\\\geq$ 2 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_2_bjets_barrel_tau"}
        )

        mutau_control_reg_2_bjets_barrel_tau = add_category(
            config,
            name="mutau_control_reg_2_bjets_barrel_tau",
            id=159 + mutau.id,
            selection=["cat_mutau"        ,
                        "ss_charge"      ,
                        "mt_cut"         ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$\\mu\\tau$ CR\n$\\\geq$ 2 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )

        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################

        mutau_signal_reg_2_bjets_no_mt = add_category(
            config,
            name="mutau_signal_reg_2_bjets_no_mt",
            id=160 + mutau.id,
            selection=["cat_mutau"        ,
                    "os_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$\\mu\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$",
            aux={'control_reg': "mutau_control_reg_2_bjets_no_mt"}
        )

        mutau_control_reg_2_bjets_no_mt = add_category(
            config,
            name="mutau_control_reg_2_bjets_no_mt",
            id=161 + mutau.id,
            selection=["cat_mutau"        ,
                    "ss_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$\\mu\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$",
        )

        mutau_signal_reg_2_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_2_bjets_endcap_tau_no_mt",
            id=162 + mutau.id,
            selection=["cat_mutau"        ,
                    "os_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$\\mu\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "mutau_control_reg_2_bjets_endcap_tau_no_mt"}
        )

        mutau_control_reg_2_bjets_endcap_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_2_bjets_endcap_tau_no_mt",
            id=163 + mutau.id,
            selection=["cat_mutau"        ,
                    "ss_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$\\mu\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )

        mutau_signal_reg_2_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_signal_reg_2_bjets_barrel_tau_no_mt",
            id=164 + mutau.id,
            selection=["cat_mutau"        ,
                        "os_charge"      ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$\\mu\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "mutau_control_reg_2_bjets_barrel_tau_no_mt"}
        )

        mutau_control_reg_2_bjets_barrel_tau_no_mt = add_category(
            config,
            name="mutau_control_reg_2_bjets_barrel_tau_no_mt",
            id=165 + mutau.id,
            selection=["cat_mutau"        ,
                        "ss_charge"      ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$\\mu\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
    ################################
    ### e-tau channel categories ###
    #################################
    elif channel=='etau':
        etau = add_category(
            config,
            name="cat_etau",
            id=3,
            selection="cat_etau",
            label="$e\\tau$ incl",
        )
        #######################################################
        ### e-tau channel categories NO b jets requirements ###
        #######################################################
        etau_signal_reg_no_mT = add_category(
            config,
            name="etau_signal_reg_no_mT",
            id=206 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$ SR\nno $m_{T}$",
            aux={'control_reg': "etau_control_reg_no_mT"}
        )
        etau_control_reg_no_mT = add_category(
            config,
            name="etau_control_reg_no_mT",
            id=207 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "deep_tau_wp" ,
                       ],
            label="$e\\tau$ CR\nno $m_{T}$",
        )
        
        etau_signal_reg = add_category(
            config,
            name="etau_signal_reg",
            id=208 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$ SR\n",
            aux={'control_reg': "etau_control_reg"}
        )
        etau_control_reg = add_category(
            config,
            name="etau_control_reg",
            id=209 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       ],
            label="$e\\tau$ CR\n",
        )
        #########################################
        ### e-tau channel categories 0 b jets ###
        #########################################
        #################################
        ### SIGNAL REGION with mT cut ###
        #################################
        
        etau_signal_reg_0_bjets = add_category(
            config,
            name="etau_signal_reg_0_bjets",
            id=166 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets",
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$ SR\n0 bjets",
            aux={'control_reg': "etau_control_reg_0_bjets"}
        )
        etau_control_reg_0_bjets = add_category(
            config,
            name="etau_control_reg_0_bjets",
            id=167 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       "Zero_b_jets" ,
                       ],
            label="$e\\tau$ CR\n0 bjets",
        )
        etau_signal_reg_0_bjets_endcap_tau = add_category(
            config,
            name="etau_signal_reg_0_bjets_endcap_tau",
            id=168 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       "OC_lepton_veto",
                       ],
            label = "$e\\tau$ SR\n0 bjets\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_0_bjets_endcap_tau"}
        )
        etau_control_reg_0_bjets_endcap_tau = add_category(
            config,
            name="etau_control_reg_0_bjets_endcap_tau",
            id=169 + etau.id,
            selection=["cat_etau"   ,
                       "ss_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
            label = "$e\\tau$ SR\n0 bjets\n$\\eta_{\\tau} > 1.2$",

        )
        
        etau_signal_reg_0_bjets_barrel_tau = add_category(
            config,
            name="etau_signal_reg_0_bjets_barrel_tau",
            id=170 + etau.id,
            selection=["cat_etau"    ,
                        "os_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n0 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_0_bjets_barrel_tau"}
        )
        
        etau_control_reg_0_bjets_barrel_tau = add_category(
            config,
            name="etau_control_reg_0_bjets_barrel_tau",
            id=171 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$e\\tau$ CR\n0 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )
        
        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################
        
        etau_signal_reg_0_bjets_no_mt = add_category(
            config,
            name="etau_signal_reg_0_bjets_no_mt",
            id=172 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                      "Zero_b_jets" ,
                      "OC_lepton_veto",
                       ],
            label="$e\\tau$\n0 bjets\nno $m_{T}$",
            aux={'control_reg': "etau_control_reg_0_bjets_no_mt"}
        )
        
        etau_control_reg_0_bjets_no_mt = add_category(
            config,
            name="etau_control_reg_0_bjets_no_mt",
            id=173 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "deep_tau_wp" ,
                       "Zero_b_jets" ,
                       ],
            label="$e\\tau$ CR\n0 bjets\nno $m_{T}$",
        )
        etau_signal_reg_0_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_0_bjets_endcap_tau_no_mt",
            id=174 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$ SR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_0_bjets_endcap_tau_no_mt"}
        )
        etau_control_reg_0_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_control_reg_0_bjets_endcap_tau_no_mt",
            id=175 + etau.id,
            selection=["cat_etau"   ,
                       "ss_charge"  ,
                       "deep_tau_wp",
                       "Zero_b_jets" ,
                       "tau_endcap" ,
                       ],
            label="$e\\tau$ CR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )
        etau_signal_reg_0_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_0_bjets_barrel_tau_no_mt",
            id=176 + etau.id,
            selection=["cat_etau"    ,
                        "os_charge"  ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_0_bjets_barrel_tau_no_mt"}
        )
        
        etau_control_reg_0_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_control_reg_0_bjets_barrel_tau_no_mt",
            id=177 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "deep_tau_wp",
                        "Zero_b_jets",
                        "tau_barrel" ,
                        ],
            label="$e\\tau$ CR\n0 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
        
        #########################################
        ### e-tau channel categories 1 b jets ###
        #########################################

        #################################
        ### SIGNAL REGION with mT cut ###
        #################################

        etau_signal_reg_1_bjets = add_category(
            config,
            name="etau_signal_reg_1_bjets",
            id=178 + etau.id,
            selection=["cat_etau"   ,
                    "os_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets",
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n1 bjet",
            aux={'control_reg': "etau_control_reg_1_bjets"}
        )

        etau_control_reg_1_bjets = add_category(
            config,
            name="etau_control_reg_1_bjets",
            id=179 + etau.id,
            selection=["cat_etau"    ,
                    "ss_charge"   ,
                    "mt_cut"      ,
                    "deep_tau_wp" ,
                    "One_b_jets"  ,
                    ],
            label="$e\\tau$ CR\n1 bjet",
        )

        etau_signal_reg_1_bjets_endcap_tau = add_category(
            config,
            name="etau_signal_reg_1_bjets_endcap_tau",
            id=180 + etau.id,
            selection=["cat_etau"   ,
                    "os_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n1 bjets\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_1_bjets_endcap_tau"}
        )

        etau_control_reg_1_bjets_endcap_tau = add_category(
            config,
            name="etau_control_reg_1_bjets_endcap_tau",
            id=181 + etau.id,
            selection=["cat_etau"   ,
                    "ss_charge"  ,
                    "mt_cut"     ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$e\\tau$ CR\n1 bjets\n$\\eta_{\\tau} > 1.2$",
        )

        etau_signal_reg_1_bjets_barrel_tau = add_category(
            config,
            name="etau_signal_reg_1_bjets_barrel_tau",
            id=182 + etau.id,
            selection=["cat_etau"    ,
                        "os_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "One_b_jets" ,
                        "tau_barrel" ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n1 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_1_bjets_barrel_tau"}
        )

        etau_control_reg_1_bjets_barrel_tau = add_category(
            config,
            name="etau_control_reg_1_bjets_barrel_tau",
            id=183 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "mt_cut"     ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$e\\tau$ CR\n1 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )

        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################

        etau_signal_reg_1_bjets_no_mt = add_category(
            config,
            name="etau_signal_reg_1_bjets_no_mt",
            id=184 + etau.id,
            selection=["cat_etau"   ,
                    "os_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n1 bjets\nno $m_{T}$",
            aux={'control_reg': "etau_control_reg_1_bjets_no_mt"}
        )

        etau_control_reg_1_bjets_no_mt = add_category(
            config,
            name="etau_control_reg_1_bjets_no_mt",
            id=185 + etau.id,
            selection=["cat_etau"    ,
                    "ss_charge"   ,
                    "deep_tau_wp" ,
                    "One_b_jets"  ,
                    ],
            label="$e\\tau$ CR\n1 bjets\nno $m_{T}$",
        )

        etau_signal_reg_1_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_1_bjets_endcap_tau_no_mt",
            id=186 + etau.id,
            selection=["cat_etau"   ,
                    "os_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_1_bjets_endcap_tau_no_mt"}
        )

        etau_control_reg_1_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_control_reg_1_bjets_endcap_tau_no_mt",
            id=187 + etau.id,
            selection=["cat_etau"   ,
                    "ss_charge"  ,
                    "deep_tau_wp",
                    "One_b_jets" ,
                    "tau_endcap" ,
                    ],
            label="$e\\tau$ CR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )

        etau_signal_reg_1_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_1_bjets_barrel_tau_no_mt",
            id=188 + etau.id,
            selection=["cat_etau"    ,
                        "os_charge"  ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_1_bjets_barrel_tau_no_mt"}
        )

        etau_control_reg_1_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_control_reg_1_bjets_barrel_tau_no_mt",
            id=189 + etau.id,
            selection=["cat_etau"  ,
                        "ss_charge"  ,
                        "deep_tau_wp",
                        "One_b_jets",
                        "tau_barrel" ,
                        ],
            label="$e\\tau$ CR\n1 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
            
        #########################################
        ### e-tau channel categories >= 2 bjets ###
        #########################################

        #################################
        ### SIGNAL REGION with mT cut ###
        #################################

        etau_signal_reg_2_bjets = add_category(
            config,
            name="etau_signal_reg_2_bjets",
            id=190 + etau.id,
            selection=["cat_etau"        ,
                    "os_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets",
            aux={'control_reg': "etau_control_reg_2_bjets"}
        )

        etau_control_reg_2_bjets = add_category(
            config,
            name="etau_control_reg_2_bjets",
            id=191 + etau.id,
            selection=["cat_etau"        ,
                    "ss_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets",
        )

        etau_signal_reg_2_bjets_endcap_tau = add_category(
            config,
            name="etau_signal_reg_2_bjets_endcap_tau",
            id=192 + etau.id,
            selection=["cat_etau"        ,
                    "os_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets\n$\\eta_{\\tau}  > 1.2$",
            aux={'control_reg': "etau_control_reg_2_bjets_endcap_tau"}
        )

        etau_control_reg_2_bjets_endcap_tau = add_category(
            config,
            name="etau_control_reg_2_bjets_endcap_tau",
            id=193 + etau.id,
            selection=["cat_etau"        ,
                    "ss_charge"       ,
                    "mt_cut"          ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets\n$\\eta_{\\tau}  > 1.2$",
        )

        etau_signal_reg_2_bjets_barrel_tau = add_category(
            config,
            name="etau_signal_reg_2_bjets_barrel_tau",
            id=194 + etau.id,
            selection=["cat_etau"        ,
                        "os_charge"      ,
                        "mt_cut"         ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_2_bjets_barrel_tau"}
        )

        etau_control_reg_2_bjets_barrel_tau = add_category(
            config,
            name="etau_control_reg_2_bjets_barrel_tau",
            id=195 + etau.id,
            selection=["cat_etau"        ,
                        "ss_charge"      ,
                        "mt_cut"         ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets\n$\\eta_{\\tau} \\leq 1.2$",
        )

        ####################################
        ### SIGNAL REGION without mT cut ###
        ####################################

        etau_signal_reg_2_bjets_no_mt = add_category(
            config,
            name="etau_signal_reg_2_bjets_no_mt",
            id=196 + etau.id,
            selection=["cat_etau"        ,
                    "os_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$",
            aux={'control_reg': "etau_control_reg_2_bjets_no_mt"}
        )

        etau_control_reg_2_bjets_no_mt = add_category(
            config,
            name="etau_control_reg_2_bjets_no_mt",
            id=197 + etau.id,
            selection=["cat_etau"        ,
                    "ss_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$",
        )

        etau_signal_reg_2_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_2_bjets_endcap_tau_no_mt",
            id=198 + etau.id,
            selection=["cat_etau"        ,
                    "os_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    "OC_lepton_veto",
                    ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
            aux={'control_reg': "etau_control_reg_2_bjets_endcap_tau_no_mt"}
        )

        etau_control_reg_2_bjets_endcap_tau_no_mt = add_category(
            config,
            name="etau_control_reg_2_bjets_endcap_tau_no_mt",
            id=199 + etau.id,
            selection=["cat_etau"        ,
                    "ss_charge"       ,
                    "deep_tau_wp"     ,
                    "At_least_2_b_jets",
                    "tau_endcap"      ,
                    ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} > 1.2$",
        )

        etau_signal_reg_2_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_signal_reg_2_bjets_barrel_tau_no_mt",
            id=200 + etau.id,
            selection=["cat_etau"        ,
                        "os_charge"      ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        "OC_lepton_veto",
                        ],
            label="$e\\tau$ SR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
            aux={'control_reg': "etau_control_reg_2_bjets_barrel_tau_no_mt"}
        )

        etau_control_reg_2_bjets_barrel_tau_no_mt = add_category(
            config,
            name="etau_control_reg_2_bjets_barrel_tau_no_mt",
            id=201 + etau.id,
            selection=["cat_etau"        ,
                        "ss_charge"      ,
                        "deep_tau_wp"    ,
                        "At_least_2_b_jets",
                        "tau_barrel"     ,
                        ],
            label="$e\\tau$ CR\n$\\geq$ 2 bjets\nno $m_{T}$\n$\\eta_{\\tau} \\leq 1.2$",
        )
        
        ##### Categories to plot Jet Variables with mT cut###
        etau_signal_reg_b_jets = add_category(
            config,
            name="etau_signal_reg_b_jets",
            id=202 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "mt_cut"     ,
                       "deep_tau_wp",
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$\nbjets incl",
            aux={'control_reg': "etau_control_reg_b_jets"}
        )
        etau_control_reg_b_jets = add_category(
            config,
            name="etau_control_reg_b_jets",
            id=203 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "mt_cut"      ,
                       "deep_tau_wp" ,
                       ],
            label="$e\\tau$ CR\nbjets",
        )
        ##### Categories to plot Jet Variables without mT cut###
        etau_signal_reg_b_jets_no_mT = add_category(
            config,
            name="etau_signal_reg_b_jets_no_mT",
            id=204 + etau.id,
            selection=["cat_etau"   ,
                       "os_charge"  ,
                       "deep_tau_wp",
                       "OC_lepton_veto",
                       ],
            label="$e\\tau$\nbjets incl\nno $m_{T}$",
            aux={'control_reg': "etau_control_reg_b_jets_no_mT"}
        )
        etau_control_reg_b_jets_no_mT = add_category(
            config,
            name="etau_control_reg_b_jets_no_mT",
            id=205 + etau.id,
            selection=["cat_etau"    ,
                       "ss_charge"   ,
                       "deep_tau_wp" ,
                       ],
            label="$e\\tau$ CR\nbjets\nno $m_{T}$",
        )
    ################################
    ### e-mu channel categories ###
    #################################
    elif channel=='emu':
        emu = add_category(
            config,
            name="cat_emu",
            id=4,
            selection="cat_emu",
            label="$e\\mu$ incl",
        )
        #######################################################
        ### e-tau channel categories NO b jets requirements ###
        #######################################################
        emu_signal_reg = add_category(
            config,
            name="emu_signal_reg",
            id=207 + emu.id,
            selection=["cat_emu"  ,
                       "os_charge",
                       ],
            label="$e\\mu$ SR\n",
            aux={'control_reg': "emu_control_reg"}
        )
        emu_control_reg = add_category(
            config,
            name="emu_control_reg",
            id=208 + emu.id,
            selection=["cat_emu"  ,
                       "ss_charge",
                       ],
            label="$e\\mu$ CR\n",
        )
        # emu_signal_reg_no_os = add_category(
        #     config,
        #     name="emu_signal_reg_no_os",
        #     id=209 + emu.id,
        #     selection=["cat_emu"  ,
        #                ],
        #     label="$e\\mu$ SR\n",
        # )
    
    elif channel=='tautau':
        tautau = add_category(
            config,
            name="cat_tautau",
            id=5,
            selection="cat_tautau",
            label="$\tau\tau$ incl",
        )

    else:
        raise ValueError(f"attempt to process more than one channel: {ch_str}")


    

  
