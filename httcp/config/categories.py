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
    add_category(
        config,
        name="2j",
        id=100,
        selection="cat_2j",
        label="2 jets",
    )

    # ------------------------------- #
    #              e-tau              #
    # ------------------------------- #
    add_category(
        config,
        name="etau",
        id=101,
        selection="sel_etau",
        label="etau_channel",
    )
    add_category(
        config,
        name="etau_pion",
        id=102,
        selection="sel_etau_pion",
        label="etau_channel_pi",
    )
    add_category(
        config,
        name="etau_rho",
        id=103,
        selection="sel_etau_rho",
        label="etau_channel_rho",
    )
    add_category(
        config,
        name="etau_a1",
        id=104,
        selection="sel_etau_a1",
        label="etau_channel_a1",
    )

    # ------------------------------- #
    #              mu-tau             #
    # ------------------------------- #
    add_category(
        config,
        name="mutau",
        id=201,
        selection="sel_mutau",
        label="mutau_channel",
    )
    add_category(
        config,
        name="mutau_pion",
        id=202,
        selection="sel_mutau_pion",
        label="mutau_channel_pi",
    )
    add_category(
        config,
        name="mutau_rho",
        id=203,
        selection="sel_mutau_rho",
        label="mutau_channel_rho",
    )
    add_category(
        config,
        name="mutau_a1",
        id=204,
        selection="sel_mutau_a1",
        label="mutau_channel_a1",
    )

    # ------------------------------- #
    #             tau-tau             #
    # ------------------------------- #
    add_category(
        config,
        name="tautau",
        id=301,
        selection="sel_tautau",
        label="tautau_channel",
    )
    add_category(
        config,
        name="tautau_pionpion",
        id=302,
        selection="sel_tautau_pionpion",
        label="tautau_channel_pi_pi",
    )
    add_category(
        config,
        name="tautau_rhorho",
        id=303,
        selection="sel_tautau_rhorho",
        label="tautau_channel_rho_rho",
    )
    add_category(
        config,
        name="tautau_a1a1",
        id=304,
        selection="sel_tautau_a1a1",
        label="tautau_channel_a1_a1",
    )
    add_category(
        config,
        name="tautau_pionrho",
        id=305,
        selection="sel_tautau_pionrho",
        label="tautau_channel_pi_rho",
    )
    add_category(
        config,
        name="tautau_a1pion",
        id=306,
        selection="sel_tautau_a1pion",
        label="tautau_channel_a1_pion",
    )
    add_category(
        config,
        name="tautau_a1rho",
        id=307,
        selection="sel_tautau_a1rho",
        label="tautau_channel_a1_rho",
    )
