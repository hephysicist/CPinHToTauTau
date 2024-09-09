# coding: utf-8

"""
Definition of categories.
"""

import order as od

from columnflow.config_util import add_category


def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_category(
        config,
        name="incl",
        id=1,
        selection="cat_incl",
        label="inclusive",
    )
    #Main categories for the three channels
    mutau = add_category(
        config,
        name="mutau",
        id=2,
        selection="mutau",
        label=r"$\mu\tau$ inclusive",
    )
    
    etau = add_category(
        config,
        name="etau",
        id=3,
        selection="etau",
        label=r"$e\tau$ inclusive",
    )
    
    tautau = add_category(
        config,
        name="tautau",
        id=4,
        selection="tautau",
        label=r"$\tau\tau$ inclusive",
    )
    os_charge = add_category(
        config,
        name="os_charge",
        id=5,
        selection="os_charge",
        label=r"$q_1\cdot q_2 < 0$",
    )
    
    ss_charge = add_category(
        config,
        name="ss_charge",
        id=6,
        selection="ss_charge",
        label=r"$q_1\cdot q_2 > 0$",
    )
    
    mT_cut = add_category(
        config,
        name="mt_cut",
        id=7,
        selection="mt_cut",
        label=r"$mT < 50$ GeV",
    )
    
    mT_inv_cut = add_category(
        config,
        name="mt_inv_cut",
        id=8,
        selection="mt_inv_cut",
        label=r"$mT \geq 50$ GeV",
    )
    

    # mu-tau categories
    mutau_signal_reg = mutau.add_category(
        name="mutau_signal_reg",
        id=10 + mutau.id,
        selection=[mutau.selection,
                   os_charge.selection,
                   mT_cut.selection],
        label=r"$\mu\tau$ signal region",
    )
    
    mutau_signal_reg_inv_mt = mutau.add_category(
        name="mutau_signal_reg_inv_mt",
        id=20 + mutau.id,
        selection=[mutau.selection,
                   os_charge.selection,
                   mT_inv_cut.selection],
        label=(r"$\mu\tau$ signal region"
               r"$m_T > 50$ GeV"),
    )
    mutau_control_reg = mutau.add_category(
        name="mutau_control_reg",
        id=30 + mutau.id,
        selection=[mutau.selection,
                   ss_charge.selection,
                   mT_cut.selection],
        label=r"$\mu\tau$ control region",
    )
    mutau_control_reg_inv_mt = mutau.add_category(
        name="mutau_control_reg_inv_mt",
        id=40 + mutau.id,
         selection=[mutau.selection,
                   ss_charge.selection,
                   mT_inv_cut.selection],
        label=(r"$\mu\tau$ control region, $m_T > 50$ GeV"),
    )
    
    mutau_signal_reg_no_mt = mutau.add_category(
        name="mutau_signal_reg_no_mt",
        id=50 + mutau.id,
        selection=[mutau.selection,
                   os_charge.selection],
        label=r"$\mu\tau$ signal region, no mT cut",
    )
    
    # e-tau categories
    etau_signal_reg = etau.add_category(
        name="etau_signal_reg",
        id=10 + etau.id,
        selection=[etau.selection,
                   os_charge.selection,
                   mT_cut.selection],
        label=r"$e\tau$ signal region",
    )
    
    etau_signal_reg_inv_mt = etau.add_category(
        name="etau_signal_reg_inv_mt",
        id=20 + etau.id,
        selection=[etau.selection,
                   os_charge.selection,
                   mT_inv_cut.selection],
        label=(r"$e\tau$ signal region"
               r"$m_T > 50$ GeV"),
    )
    etau_control_reg = etau.add_category(
        name="etau_control_reg",
        id=30 + etau.id,
        selection=[etau.selection,
                   ss_charge.selection,
                   mT_cut.selection],
        label=r"$e\tau$ control region",
    )
    etau_control_reg_inv_mt = etau.add_category(
        name="etau_control_reg_inv_mt",
        id=40 + etau.id,
         selection=[etau.selection,
                   ss_charge.selection,
                   mT_inv_cut.selection],
        label= (r"$e\tau$ control region"
                r"$m_T > 50$ GeV"),
    )
    
    # tau-tau categories
    tautau_signal_reg = tautau.add_category(
        name="tautau_signal_reg",
        id=10 + tautau.id,
        selection=[tautau.selection,
                   os_charge.selection,
                   mT_cut.selection],
        label=r"$\tau\tau$ signal region",
    )
        
    
    

    
    
    
    
    
    
    # add_category(
    #     config,
    #     name="2j",
    #     id=100,
    #     selection="cat_2j",
    #     label="2 jets",
    # )

    # # ------------------------------- #
    # #              e-tau              #
    # # ------------------------------- #
    # add_category(
    #     config,
    #     name="etau",
    #     id=101,
    #     selection="sel_etau",
    #     label="etau_channel",
    # )
    # add_category(
    #     config,
    #     name="etau_pion",
    #     id=102,
    #     selection="sel_etau_pion",
    #     label="etau_channel_pi",
    # )
    # add_category(
    #     config,
    #     name="etau_rho",
    #     id=103,
    #     selection="sel_etau_rho",
    #     label="etau_channel_rho",
    # )
    # add_category(
    #     config,
    #     name="etau_a1",
    #     id=104,
    #     selection="sel_etau_a1",
    #     label="etau_channel_a1",
    # )

    # # ------------------------------- #
    # #              mu-tau             #
    # # ------------------------------- #
    # add_category(
    #     config,
    #     name="mutau",
    #     id=201,
    #     selection="sel_mutau",
    #     label="mutau_channel",
    # )
    # add_category(
    #     config,
    #     name="mutau_pion",
    #     id=202,
    #     selection="sel_mutau_pion",
    #     label="mutau_channel_pi",
    # )
    # add_category(
    #     config,
    #     name="mutau_rho",
    #     id=203,
    #     selection="sel_mutau_rho",
    #     label="mutau_channel_rho",
    # )
    # add_category(
    #     config,
    #     name="mutau_a1",
    #     id=204,
    #     selection="sel_mutau_a1",
    #     label="mutau_channel_a1",
    # )

    # # ------------------------------- #
    # #             tau-tau             #
    # # ------------------------------- #
    # add_category(
    #     config,
    #     name="tautau",
    #     id=301,
    #     selection="sel_tautau",
    #     label="tautau_channel",
    # )
    # add_category(
    #     config,
    #     name="tautau_pionpion",
    #     id=302,
    #     selection="sel_tautau_pionpion",
    #     label="tautau_channel_pi_pi",
    # )
    # add_category(
    #     config,
    #     name="tautau_rhorho",
    #     id=303,
    #     selection="sel_tautau_rhorho",
    #     label="tautau_channel_rho_rho",
    # )
    # add_category(
    #     config,
    #     name="tautau_a1a1",
    #     id=304,
    #     selection="sel_tautau_a1a1",
    #     label="tautau_channel_a1_a1",
    # )
    # add_category(
    #     config,
    #     name="tautau_pionrho",
    #     id=305,
    #     selection="sel_tautau_pionrho",
    #     label="tautau_channel_pi_rho",
    # )
    # add_category(
    #     config,
    #     name="tautau_a1pion",
    #     id=306,
    #     selection="sel_tautau_a1pion",
    #     label="tautau_channel_a1_pion",
    # )
    # add_category(
    #     config,
    #     name="tautau_a1rho",
    #     id=307,
    #     selection="sel_tautau_a1rho",
    #     label="tautau_channel_a1_rho",
    # )
