# coding: utf-8

"""
Style definitions.
"""

from __future__ import annotations

import re
import random
from collections import defaultdict

import order as od

from columnflow.util import DotDict, try_int
from columnflow.types import Callable

def setup_plot_styles(config: od.Config) -> None:
    """
    Setup plot styles.
    """
    # general settings
    config.x.default_general_settings = {
        "cms_label": "wip", "whitespace_fraction": 0.31,
    }

    # default component configs
    gridspec = {
        "height_ratios": [3, 0.9],
    }
    legend = {
        "borderpad": 0, "borderaxespad": 1.2, "columnspacing": 1.8, "labelspacing": 0.28, "fontsize": 16,
        "cf_line_breaks": True, "cf_short_labels": False,
    }
    ratio = {
        "yloc": "center",
    }
    annotate = {
        "fontsize": 18, "style": "italic", "xycoords": "axes fraction", "xy": (0.035, 0.955),
    }

    # wide legend
    # - 3 columns, backgrounds in first 2 columns
    # - shortened process labels
    # - changed annotation (channel) position to fit right under legend
    wide_legend = legend | {
        "ncols": 2, "loc": "upper right", "cf_entries_per_column": legend_entries_per_column, "cf_short_labels": True,
    }
    annotate_wide = annotate | {
        "xy": (0.035, 0.765),
    }

    # wide extended legend, same as wide legend except
    # - process labels are not shortened
    # - annotation (channel) moved slightly down to fut under (now taller) legend
    wide_ext_legend = wide_legend | {
        "cf_short_labels": False,
    }
    annotate_wide_ext = annotate_wide | {
        "xy": (0.035, 0.750),
    }

    # construct named style configs
    config.x.custom_style_config_groups = {
        "default": (default_cfg := {
            "gridspec_cfg": gridspec,
            "rax_cfg": ratio,
            "legend_cfg": legend,
            "annotate_cfg": annotate,
        }),
        "wide_legend": (wide_legend_cfg := {
            **default_cfg,
            "legend_cfg": wide_legend,
            "annotate_cfg": annotate_wide,
        }),
        "wide_ext_legend": {
            **wide_legend_cfg,
            "legend_cfg": wide_ext_legend,
            "annotate_cfg": annotate_wide_ext,
        },
    }

    config.x.default_custom_style_config = "wide_legend"
    #config.x.default_blinding_threshold = 0


def legend_entries_per_column(ax, handles: list, labels: list, n_cols: int) -> list[int]:
    """
    Controls number of entries such that backgrounds are in the first n - 1 columns, and everything
    else in the last one.
    """
    # get number of background and remaining entries
    n_backgrounds = sum(1 for handle in handles if handle.__class__.__name__ == "StepPatch")
    n_other = len(handles) - n_backgrounds

    # fill number of entries per column
    entries_per_col = n_cols * [0]
    n_bkg_cols = n_cols
    # set last column if non-backgrounds are present
    if n_other:
        entries_per_col[-1] = n_other
        n_bkg_cols -= 1
    # fill background columns
    for i in range(n_bkg_cols):
        entries_per_col[i] = n_backgrounds // n_bkg_cols + (n_backgrounds % n_bkg_cols > i)

    return entries_per_col
