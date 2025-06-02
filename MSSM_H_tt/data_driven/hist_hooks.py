"""
Histogram hooks.
"""

from collections import defaultdict

import law
import order as od
import scinum as sn

from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
hist = maybe_import("hist")


logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    # helper to convert a histogram to a number object containing bin values and uncertainties
    # from variances stored in an array of values
    def hist_to_num(h: hist.Hist, unc_name=str(sn.DEFAULT)) -> sn.Number:
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})

    # # helper to integrate values stored in an array based number object
    # def integrate_num(num: sn.Number, axis=None) -> sn.Number:
    #     return sn.Number(
    #         nominal=num.nominal.sum(axis=axis),
    #         uncertainties={
    #             unc_name: (
    #                 (unc_values_up**2).sum(axis=axis)**0.5,
    #                 (unc_values_down**2).sum(axis=axis)**0.5,
    #             )
    #             for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
    #         },
    #     )

    def qcd_estimation(task, hists):
        if not hists:
            return hists
        from cmsdb.processes.qcd import qcd
        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups
        signal_categories : dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        #Filling the dictionay with signal histograms
        for the_name in config.categories.names():
            if 'signal_reg' in the_name:
                signal_categories[the_name] = config.get_category(the_name)
        

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd] = qcd_hist = mc_hist.copy().reset()
        
        def get_hist (h, category): 
            return h[{"category": hist.loc(category.id)}]
        
        
        def find_by_id(hist, cat_id):
            cat_axis = hist.axes["category"]
            for  the_idx, the_id in enumerate(cat_axis): 
                if the_id == cat_id: return the_idx
            return -1

        for (the_name, the_category) in signal_categories.items():
            control_cat = config.get_category(the_category.aux['control_reg'])
            data_array = hist_to_num(get_hist(data_hist, control_cat))
            mc_array = hist_to_num(get_hist(mc_hist, control_cat))
            qcd_array = data_array - mc_array
            qcd_values = qcd_array()
            qcd_values[qcd_array() < 0] = 0
            qcd_variances = qcd_array(sn.UP, sn.ALL, unc=True)**2     
            cat_idx = find_by_id(qcd_hist, the_category.id)
            QCD_weight = 1
            qcd_hist.view().value[cat_idx, ...] = qcd_values*QCD_weight
            qcd_hist.view().variance[cat_idx, ...] = qcd_variances*QCD_weight
          
        # for group_name in complete_groups:
        
        #     group = qcd_groups[group_name]

        #     # get the corresponding histograms and convert them to number objects,
        #     # each one storing an array of values with uncertainties
        #     # shapes: (SHIFT, VAR)
        #     get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
        #     os_noniso_mc = hist_to_num(get_hist(mc_hist, "os_noniso"), "os_noniso_mc")
        #     ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
        #     ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
        #     os_noniso_data = hist_to_num(get_hist(data_hist, "os_noniso"), "os_noniso_data")
        #     ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
        #     ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

        #     # estimate qcd shapes in the three sideband regions
        #     # shapes: (SHIFT, VAR)
        #     os_noniso_qcd = os_noniso_data - os_noniso_mc
        #     ss_iso_qcd = ss_iso_data - ss_iso_mc
        #     ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

        #     # get integrals in ss regions for the transfer factor
        #     # shapes: (SHIFT,)
        #     int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
        #     int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)

        #     # complain about negative integrals
        #     int_ss_iso_neg = int_ss_iso <= 0
        #     int_ss_noniso_neg = int_ss_noniso <= 0
        #     if int_ss_iso_neg.any():
        #         shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_iso_neg)[0]))
        #         shifts = list(map(config.get_shift, shift_ids))
        #         logger.warning(
        #             f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
        #             f"{', '.join(map(str, shifts))}",
        #         )
        #     if int_ss_noniso_neg.any():
        #         shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_noniso_neg)[0]))
        #         shifts = list(map(config.get_shift, shift_ids))
        #         logger.warning(
        #             f"negative QCD integral in ss_noniso region for group {group_name} and shifts: "
        #             f"{', '.join(map(str, shifts))}",
        #         )

        #     # ABCD method
        #     # shape: (SHIFT, VAR)
        #     os_iso_qcd = os_noniso_qcd * ((int_ss_iso / int_ss_noniso)[:, None])

        #     # combine uncertainties and store values in bare arrays
        #     os_iso_qcd_values = os_iso_qcd()
        #     os_iso_qcd_variances = os_iso_qcd(sn.UP, sn.ALL, unc=True)**2

        #     # define uncertainties
        #     unc_data = os_iso_qcd(sn.UP, ["os_noniso_data", "ss_iso_data", "ss_noniso_data"], unc=True)
        #     unc_mc = os_iso_qcd(sn.UP, ["os_noniso_mc", "ss_iso_mc", "ss_noniso_mc"], unc=True)
        #     unc_data_rel = abs(unc_data / os_iso_qcd_values)
        #     unc_mc_rel = abs(unc_mc / os_iso_qcd_values)

        #     # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
        #     keep_variance_mask = (
        #         np.isfinite(unc_mc_rel) &
        #         (unc_mc_rel > unc_data_rel) &
        #         (unc_mc_rel > 0.15)
        #     )
        #     os_iso_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
        #     os_iso_qcd_variances[~keep_variance_mask] = 0

        #     # retro-actively set values to zero for shifts that had negative integrals
        #     neg_int_mask = int_ss_iso_neg | int_ss_noniso_neg
        #     os_iso_qcd_values[neg_int_mask] = 1e-5
        #     os_iso_qcd_variances[neg_int_mask] = 0

        #     # residual zero filling
        #     zero_mask = os_iso_qcd_values <= 0
        #     os_iso_qcd_values[zero_mask] = 1e-5
        #     os_iso_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
        return hists

    config.x.hist_hooks = {
        "qcd_hook": qcd_estimation,
    }
