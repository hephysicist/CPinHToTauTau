# coding: utf-8

"""
Histogram hooks.
"""

import law
import order as od
import scinum as sn

from columnflow.util import maybe_import
from law.util import InsertableDict

np = maybe_import("numpy")
ak = maybe_import("awkward")
hist = maybe_import("hist")

logger = law.logger.get_logger(__name__)


def calc_yields(hists: dict, locator: dict, fake_locator: dict) -> dict:
    data_hists = [h for p, h in hists.items() if p.is_data]
    data_hist = sum(data_hists[1:], data_hists[0].copy())
    data = data_hist[locator].values()

    mc_hists = [h for p, h in hists.items() if (p.is_mc and not p.has_tag("signal"))]
    mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
    mc = mc_hist[locator].values()
    mc_fakes = mc_hist[fake_locator].values()

    wj_ratio = mc_fakes / np.maximum(data - mc, 1)
    wj_err = np.ones_like(wj_ratio)

    qcd_ratio = np.maximum((data - (mc_fakes + mc)), 0) / np.maximum(data - mc, 1)
    qcd_err = np.ones_like(wj_ratio)

    return {
        "wj": wj_ratio,
        "wj_err": wj_err,
        "qcd": qcd_ratio,
        "qcd_err": qcd_err,
    }


def get_data_hist(hists: dict):
    data_hists = [h for p, h in hists.items() if p.is_data]
    if len(data_hists):
        return sum(data_hists[1:], data_hists[0].copy())
    return None


def get_mc_hist(hists: dict, remove_wj: bool = False):
    mc_hists = [
        h
        for p, h in hists.items()
        if p.is_mc
        and not p.has_tag("signal")
        and ((remove_wj and ("wj" not in p.name)) or not remove_wj)
    ]
    if len(mc_hists):
        return sum(mc_hists[1:], mc_hists[0].copy())
    return None


def get_signal_hists(hists: dict):
    return [
        h
        for p, h in hists.items()
        if p.is_mc and p.has_tag("signal") and (p.name != "qcd")
    ]


def rel_err(h_arr=None, err_arr=None):
    h_arr = h_arr or []
    err_arr = err_arr or []
    if len(h_arr):
        sum_var = np.zeros_like(h_arr[0].values())
    else:
        sum_var = np.zeros_like(err_arr[0])
    for x in h_arr:
        sum_var += x.variances() / np.maximum(x.values() ** 2, 1)
    for the_arr in err_arr:
        sum_var += the_arr
    return np.sqrt(sum_var)


def find_idxs(h: hist.Hist, cat: str, shift: str):
    return h.axes[:-1].index(cat, shift)


def drop_categories(h, axis_name, keep):
    axes_new = [
        hist.axis.StrCategory(keep, name=ax.name) if ax.name == axis_name else ax
        for ax in h.axes
    ]

    h_new = hist.Hist(*axes_new)

    for cat in keep:
        h_new[{axis_name: cat}] = h[{axis_name: cat}]

    return h_new


def add_hist_hooks(analysis: "od.Analysis") -> None:
    """
    Add histogram hooks to a configuration.
    """
   
    flat_tf = True
    # def qcd_estimation(task, inputs, variable_name=None, category_name=None, **kwargs): #cf0p3
    #     output = {}
    #     for config, hists in inputs.items():
    #         if task.get_task_family() == 'cf.CreateDatacards':
    #             sr_cats = [task.branch_data.config_data[config.name].category]
    #             incl_h = sum(list(hists.values()))
    #             ax = incl_h.axes['shift']
    #             shifts = [ax.value(i) for i in range(ax.size)]
    #         else: 
    #             sr_cats = []
    #             if 'incl' in task.hist_hooks:
    #                 print("Updating category list to create cats inclusive by tau decay channels")
    #                 decay_ch = ['tau2pi','tau2rho','tau2a1','tau2a1_3pr']
    #                 for the_cat in task.categories:
    #                     for the_ch in decay_ch:
    #                         if '_tau2' not in the_cat:
    #                             sr_cats.append( '__'.join((the_cat,the_ch)))
    #                         else:
    #                             sr_cats.append(the_cat)
    #             else:
    #                 sr_cats = task.categories
    #             shifts = getattr(task, "shifts", ["nominal"])
    #         for shift in shifts:
    #             for the_cat in sr_cats:
    #                 print(f'producing qcd for {the_cat}, shift: {shift}')
    #                 sr = config.get_category(the_cat)
    #                 d = {}
    #                 mc = {}
    #                 cr_cat = ''
    #                 for reg_name, full_name in sr.aux['abcd_regs'].items():   
    #                     print("WARNING: using nominal shift for qcd estimation")
    #                     loc_dict = {'category': hist.loc(full_name),'shift': hist.loc('nominal')}
    #                     full_d = get_data_hist(hists)
    #                     full_mc = get_mc_hist(hists)   
    #                     if full_name in list(full_d.axes[0]): d[reg_name] = get_data_hist(hists)[loc_dict]
    #                     if full_name in list(full_mc.axes[0]): mc[reg_name] = get_mc_hist(hists)[loc_dict]
    #                     if reg_name == 'dr_num': cr_cat = full_name
                            
    #                 from cmsdb.processes.qcd import qcd
    #                 h_donor_name  = list(hists.keys())[0]
    #                 if qcd not in hists.keys():
    #                     hists[qcd] = hists[h_donor_name].copy().reset()
    #                     tmp_arr = hists[qcd].view()
    #                 if 'ar' in d.keys() and not d['ar'].empty():
    #                     if flat_tf:
    #                         tf = 1
    #                         #num = ak.sum(d['dr_num'].values() - mc['dr_num'].values())
    #                         #den = ak.sum(d['dr_den'].values() - mc['dr_den'].values())
    #                         #if (num > 0) and (den > 0):
    #                         #    tf = num/den
    #                         #else:
    #                         #    tf = 1. 
                        
    #                     else:
    #                         num = d['dr_num'].values() - mc['dr_num'].values()
    #                         den = d['dr_den'].values() - mc['dr_den'].values()

    # flat_tf = True
    
    
    
    #          h_qcd[...] = tmp_arr
    #          output[config] = hists

    #      return output

    def _get_config_process(config, process_name: str, required: bool = False):
        try:
            return config.get_process(process_name)
        except Exception:
            if required:
                raise
            return None

    def _get_or_create_hist_like(hists: dict, proc_obj, donor_hist: hist.Hist, shift_sources):
        if proc_obj in hists:
            return hists[proc_obj]

        def storage_from(donor: hist.Hist):
            for attr in ("_storage_type", "storage_type"):
                if hasattr(donor, attr):
                    st = getattr(donor, attr)
                    try:
                        st = st() if callable(st) else st
                    except TypeError:
                        pass
                    try:
                        return st() if isinstance(st, type) else st
                    except Exception:
                        return st
            return hist.storage.Weight()

        axes = []
        for ax in donor_hist.axes:
            if getattr(ax, "name", None) == "shift":
                axes.append(hist.axis.StrCategory(list(shift_sources), name="shift", growth=True))
            else:
                axes.append(ax)

        hists[proc_obj] = hist.Hist(*axes, storage=storage_from(donor_hist)).reset()
        return hists[proc_obj]

    def qcd_estimation(task, inputs, variable_name=None, category_name=None, **kwargs):
        """
        QCD estimation hook compatible with both normal histogram tasks and CreateDatacards.
        """

        output = {}

        def up_down_once(sources):
            bases = []
            for s in sources:
                if s.endswith("_down"):
                    s = s[:-5]
                elif s.endswith("_up"):
                    s = s[:-3]
                bases.append(s)

            out = []
            seen = set()
            for b in bases:
                if b in seen:
                    continue
                seen.add(b)
                out.extend((f"{b}_down", f"{b}_up"))
            if "nominal" not in out:
                out.append("nominal")
            return tuple(out)

        def fmt_vals(a):
            a = np.asarray(a)
            return np.array2string(a, precision=6, separator=", ", threshold=a.size)

        for config, hists in inputs.items():
            if not hists:
                output[config] = hists
                continue

            qcd_proc = _get_config_process(config, "qcd")
            if qcd_proc is None:
                logger.warning(
                    f"config {config.name} has no process 'qcd'; skipping qcd hist hook"
                )
                output[config] = hists
                continue

            if task.get_task_family() == "cf.CreateDatacards":
                if category_name is None:
                    raise ValueError(
                        "qcd_estimation requires 'category_name' when run from CreateDatacards"
                    )
                sr_cats = [category_name]

                incl_h = sum(list(hists.values())[1:], list(hists.values())[0].copy())
                ax = incl_h.axes["shift"]
                shift_sources = tuple(ax.value(i) for i in range(ax.size))
            else:
                sr_cats = []
                if 'incl' in task.hist_hooks:
                    print("Updating category list to create cats inclusive by tau decay channels")
                    decay_ch = ['tau2pi','tau2rho','tau2a1','tau2a1_3pr']
                    for the_cat in task.categories:
                        for the_ch in decay_ch:
                            if '_tau2' not in the_cat:
                                sr_cats.append( '__'.join((the_cat,the_ch)))
                            else:
                                sr_cats.append(the_cat)
                else:
                    sr_cats = task.categories
                shift_sources = up_down_once(getattr(task, "shifts", ["nominal"]))

            full_d = get_data_hist(hists)
            full_mc = get_mc_hist(hists)

            if full_d is None:
                logger.warning(
                    f"no data histogram found for config {config.name}, skipping qcd estimation"
                )
                output[config] = hists
                continue

            if full_mc is None:
                logger.warning(
                    f"no MC histogram found for config {config.name}, skipping qcd estimation"
                )
                output[config] = hists
                continue

            def mc_shift(shift: str) -> str:
                try:
                    return shift if shift in list(full_mc.axes["shift"]) else "nominal"
                except Exception:
                    return "nominal"

            h_donor_proc = list(hists.keys())[0]
            h_qcd = _get_or_create_hist_like(
                hists=hists,
                proc_obj=qcd_proc,
                donor_hist=hists[h_donor_proc],
                shift_sources=shift_sources,
            )
            tmp_arr = h_qcd.view()

            for the_cat in sr_cats:
                print(f"producing qcd for {the_cat}, {config.name}")
                sr = config.get_category(the_cat)

                if "abcd_regs" not in sr.aux:
                    logger.warning(
                        f"category '{the_cat}' in config '{config.name}' has no aux['abcd_regs']; "
                        "skipping qcd estimation"
                    )
                    continue

                d = {}
                cr_cat = ""
                for reg_name, full_name in sr.aux["abcd_regs"].items():
                    loc_dict_data = {
                        "category": hist.loc(full_name),
                        "shift": hist.loc("nominal"),
                    }
                    if full_name in list(full_d.axes["category"]):
                        d[reg_name] = full_d[loc_dict_data]
                    if reg_name == "dr_num":
                        cr_cat = full_name
    
                if ("ar" not in d) or d["ar"].empty():
                    print("*** WARNING: AR data histogram doesn't exist or is empty! ***")
                    continue

                shift_mc_nom = mc_shift("nominal")

                mc_nom = {}
                for reg_name, full_name in sr.aux["abcd_regs"].items():
                    loc_dict_mc_nom = {
                        "category": hist.loc(full_name),
                        "shift": hist.loc(shift_mc_nom),
                    }
                    if full_name in list(full_mc.axes["category"]):
                        mc_nom[reg_name] = full_mc[loc_dict_mc_nom]

                if flat_tf:
                    tf_nom = 1.0
                else:
                    num_nom = d["dr_num"].values() - mc_nom["dr_num"].values()
                    den_nom = d["dr_den"].values() - mc_nom["dr_den"].values()
                    tf_nom = np.divide(
                        num_nom,
                        den_nom,
                        out=np.zeros_like(num_nom, dtype=float),
                        where=(den_nom != 0),
                    )

                d_ar = d["ar"]
                if "ar" in mc_nom:
                    mc_ar_val_nom = mc_nom["ar"].values()
                    mc_ar_var_nom = mc_nom["ar"].view().variance
                else:
                    mc_ar_val_nom = np.zeros_like(d_ar.values())
                    mc_ar_var_nom = np.zeros_like(d_ar.view().variance)

                val_nom = np.maximum(d_ar.values() - mc_ar_val_nom, 0.0) * tf_nom
                var_nom = (d_ar.view().variance + mc_ar_var_nom) * (tf_nom ** 2)

                cr_val_nom = None
                cr_var_nom = None
                if cr_cat:
                    if ("dr_den" in d) and ("dr_den" in mc_nom):
                        cr_val_nom = np.maximum(
                            d["dr_den"].values() - mc_nom["dr_den"].values(),
                            0.0,
                        )
                        cr_var_nom = (
                            d["dr_den"].view().variance + mc_nom["dr_den"].view().variance
                        )
                    elif "dr_den" in d:
                        cr_val_nom = d["dr_den"].values()
                        cr_var_nom = d["dr_den"].view().variance
                    else:
                        cr_val_nom = 0.0
                        cr_var_nom = 0.0

                idx_nom = find_idxs(h_qcd, the_cat, "nominal")
                tmp_arr[idx_nom].value = val_nom
                tmp_arr[idx_nom].variance = var_nom

                if cr_cat and cr_val_nom is not None:
                    idx_cr_nom = find_idxs(h_qcd, cr_cat, "nominal")
                    tmp_arr[idx_cr_nom].value = cr_val_nom
                    tmp_arr[idx_cr_nom].variance = cr_var_nom

                for shift in shift_sources:
                    if shift == "nominal":
                        continue

                    shift_mc = mc_shift(shift)

                    mc = {}
                    for reg_name, full_name in sr.aux["abcd_regs"].items():
                        loc_dict_mc = {
                            "category": hist.loc(full_name),
                            "shift": hist.loc(shift_mc),
                        }
                        if full_name in list(full_mc.axes["category"]):
                            mc[reg_name] = full_mc[loc_dict_mc]

                    if flat_tf:
                        tf = 1.0
                    else:
                        num = d["dr_num"].values() - mc["dr_num"].values()
                        den = d["dr_den"].values() - mc["dr_den"].values()
                        tf = np.divide(
                            num,
                            den,
                            out=np.zeros_like(num, dtype=float),
                            where=(den != 0),
                        )

                    if "ar" in mc:
                        mc_ar_val = mc["ar"].values()
                        mc_ar_var = mc["ar"].view().variance
                    else:
                        mc_ar_val = np.zeros_like(d_ar.values())
                        mc_ar_var = np.zeros_like(d_ar.view().variance)

                    val_before = np.maximum(d_ar.values() - mc_ar_val, 0.0) * tf
                    var_before = (d_ar.view().variance + mc_ar_var) * (tf ** 2)

                    if shift.endswith("_up") or shift.endswith("_down"):
                        print(f"[QCD] {config.name} {the_cat} {shift}: {fmt_vals(val_before)}")

                    idx = find_idxs(h_qcd, the_cat, shift)
                    tmp_arr[idx].value = val_before
                    tmp_arr[idx].variance = var_before

                    if cr_cat:
                        if ("dr_den" in d) and ("dr_den" in mc):
                            cr_val_before = np.maximum(
                                d["dr_den"].values() - mc["dr_den"].values(),
                                0.0,
                            )
                            cr_var_before = (
                                d["dr_den"].view().variance + mc["dr_den"].view().variance
                            )
                        elif "dr_den" in d:
                            cr_val_before = d["dr_den"].values()
                            cr_var_before = d["dr_den"].view().variance
                        else:
                            cr_val_before = 0.0
                            cr_var_before = 0.0

                        if shift.endswith("_up") or shift.endswith("_down"):
                            print(f"[QCD-CR] {config.name} {cr_cat} {shift}: {fmt_vals(cr_val_before)}")

                        idx_cr = find_idxs(h_qcd, cr_cat, shift)
                        tmp_arr[idx_cr].value = cr_val_before
                        tmp_arr[idx_cr].variance = cr_var_before

            h_qcd[...] = tmp_arr
            output[config] = hists
        return output

    def ff_method(task, inputs, variable_name=None, category_name=None, **kwargs):
        from cmsdb.processes.qcd import jet_fakes as jet_fakes_fallback
        output = {}
        for config, hists in inputs.items():
            qcd_proc = _get_config_process(config, "qcd")
            if qcd_proc is None:
                logger.warning(f"config {config.name} has no process 'qcd'; skipping ff_method")
                output[config] = hists
                continue

            jet_fakes_proc = _get_config_process(config, "jet_fakes") or jet_fakes_fallback

            if task.get_task_family() == "cf.CreateDatacards":
                if category_name is None:
                    raise ValueError("ff_method requires 'category_name' when run from CreateDatacards")
                print(f'Running estimation of jet fakes for {category_name}')
                sr_cats = [category_name]
                incl_h = sum(list(hists.values())[1:], list(hists.values())[0].copy())
                ax = incl_h.axes["shift"]
                var2dm = {
                    'phi_cp_mu_pi': 0,
                    'phi_cp_mu_rho': 1,
                    'phi_cp_mu_a1_1pr':2,
                    'phi_cp_mu_a1_3pr_pv_gef': 10,
                }
                
                if "phi_cp" in variable_name:
                    shifts = [s for s in ax if (("ff_" in s) and (f"dm_{var2dm[variable_name]}" in s))]  
                else:
                    shifts = [s for s in ax if "ff_" in s]  
                shifts += ["nominal"]   
            else:
                sr_cats = []
                if 'incl' in task.hist_hooks:
                    print("Updating category list to create cats inclusive by tau decay channels")
                    decay_ch = ['tau2pi','tau2rho','tau2a1','tau2a1_3pr']
                    for the_cat in task.categories:
                        for the_ch in decay_ch:
                            if '_tau2' not in the_cat:
                                sr_cats.append( '__'.join((the_cat,the_ch)))
                            else:
                                sr_cats.append(the_cat)
                else:
                    sr_cats = task.categories
                shifts = getattr(task, "shifts", ["nominal"])
            dy_proc = [p for p in hists.keys() if p.name == "dy_tt_m50"][0]
            if qcd_proc not in hists:
                hists[qcd_proc] = hists[dy_proc].copy().reset()
            tmp_qcd = hists[qcd_proc].view()
            if jet_fakes_proc not in hists:
                hists[jet_fakes_proc] = hists[dy_proc].copy().reset()
            tmp_fakes = hists[jet_fakes_proc].view()
            data = get_data_hist(hists)
            mc = get_mc_hist(hists)
            for the_cat in sr_cats:
                sr = config.get_category(the_cat)
                sr.label = sr.label.replace("prompt lep.", "") + "\njet fakes from FF"
                for the_shift in shifts:
                    print(f"Applying fake factor method on {config.name} {the_cat} {the_shift}")
                    def locator(reg, shift): 
                        return {"category": hist.loc(sr.aux["ff_regs"][reg]), "shift": hist.loc(shift)}
                    def fake_locator(reg, shift):
                        ar_string = sr.aux["ff_regs"][reg]
                        af_fakes_string = ar_string.replace("prompt", "jet_fakes")
                        return {"category": hist.loc(af_fakes_string), "shift": hist.loc(shift)}
                    qcd_shift = "nominal" if "ff_qcd" not in the_shift else the_shift
                    wj_shift = "nominal" if "ff_wj" not in the_shift else the_shift
                    data_qcd = data[locator("ar_qcd",qcd_shift)]
                    mc_qcd = mc[locator("ar_qcd",qcd_shift)]
                    data_wj = data[locator("ar_wj",wj_shift)]
                    mc_wj = mc[locator("ar_wj",wj_shift)]
                    yields = calc_yields(hists,
                                         locator("ar_yields","nominal"), 
                                         fake_locator("ar_yields","nominal"))
                    
                    tmp_qcd[find_idxs(hists[qcd_proc], the_cat, the_shift)].value = (
                        np.maximum(data_qcd.values() - mc_qcd.values(), 0) * yields["qcd"]
                    )
                    tmp_qcd[find_idxs(hists[qcd_proc], the_cat, the_shift)].variance = (
                        data_qcd.variances() + mc_qcd.variances()
                    )

                    tmp_fakes[find_idxs(hists[jet_fakes_proc], the_cat, the_shift)].value = (
                        np.maximum(data_wj.values() - mc_wj.values(), 0) * yields["wj"]
                    )
                    tmp_fakes[find_idxs(hists[jet_fakes_proc], the_cat, the_shift)].variance = (
                        data_wj.variances() + mc_wj.variances()
                    )
            hists[qcd_proc][...] = tmp_qcd.copy()
            
            wj_proc = [p for p in hists.keys() if p.name == "w"][0]
            if task.get_task_family() == "cf.CreateDatacards":
                del hists[jet_fakes_proc]
                del hists[wj_proc]
                hists[wj_proc] = hists[dy_proc].copy().reset()
                hists[wj_proc][...] = tmp_fakes.copy()
                jet_fakes_proc = wj_proc
                hists[jet_fakes_proc][...] = tmp_fakes.copy()
            else:
                del hists[wj_proc]  
                hists[jet_fakes_proc][...] = tmp_fakes.copy()
                print("Imputing nominal hists to ff variation of the SR:")
                data_proc = _get_config_process(config, "data")
                for the_proc, the_h in hists.items():
                    print(f"proc: {the_proc.name}")
                    if the_proc == data_proc:
                        hists[the_proc] = the_h[{'shift':[hist.loc('nominal')]}].copy()
                    else:
                        if the_proc == qcd_proc: 
                            skip_prfix = 'ff_wj'
                            print('QCD process: remove wj systematics')
                        elif the_proc == jet_fakes_proc: 
                            skip_prfix = 'ff_qcd'
                            print('Wj process: remove qcd systematics')
                        else: 
                            skip_prfix = 'ff_'
                        
                        skimmed_shift_axis = [s for s in shifts if skip_prfix not in s]
                        hists[the_proc] = the_h[{'shift': skimmed_shift_axis}].copy()
            output[config] = hists
        return output

    def ff_closure_test(task, inputs, variable_name=None, category_name=None, **kwargs):
        from cmsdb.processes.qcd import jet_fakes as jet_fakes_fallback

        output = {}
        for config, hists in inputs.items():
            qcd_proc = _get_config_process(config, "qcd")
            if qcd_proc is None:
                logger.warning(f"config {config.name} has no process 'qcd'; skipping ff_closure_test")
                output[config] = hists
                continue

            jet_fakes_proc = _get_config_process(config, "jet_fakes") or jet_fakes_fallback
            dr_num_cats = list(task.categories)
            if 'incl' in task.hist_hooks:
                
                    print("Updating category list to create cats inclusive by tau decay channels")
                    decay_ch = ['tau2pi','tau2rho','tau2a1','tau2a1_3pr']
                    for the_cat in task.categories:
                        for the_ch in decay_ch:
                            if '_tau2' not in the_cat:
                                dr_num_cats.append( '__'.join((the_cat,the_ch)))
                            else:
                                dr_num_cats.append(the_cat)
            shifts = getattr(task, "shifts",['nominal'])
            
            for the_cat in dr_num_cats:
                print(f"Performing closure test for {the_cat}")
                sr_name = the_cat.replace("dr_num_wj", "sr").replace("dr_num_qcd", "sr")
                sr = config.get_category(sr_name)
                data = get_data_hist(hists)
                mc = get_mc_hist(hists)
                locator = lambda reg: {"category": hist.loc(sr.aux["ff_regs"][reg]), "shift": hist.loc(task.shift)}
                iter_dict = {"wj": jet_fakes_proc, "qcd": qcd_proc}
                h_donor_name = list(hists.keys())[0]

                for name, proc in iter_dict.items():
                    if proc not in hists:
                        hists[proc] = hists[h_donor_name].copy().reset()

                    tmp_fakes = hists[proc].view()
                    h_d = data[locator(f"dr_den_{name}_w_ff")]
                    h_mc = mc[locator(f"dr_den_{name}_w_ff")]
                    num = sr.aux["ff_regs"][f"dr_num_{name}"]

                    tmp_fakes[find_idxs(hists[proc], num, task.shift)].value = np.maximum(
                        h_d.values() - h_mc.values(),
                        0,
                    )
                    tmp_fakes[find_idxs(hists[proc], num, task.shift)].variance = (
                        h_d.variances() + h_mc.variances()
                    )
                    if name =='wj':
                        wj_proc = [p for p in hists.keys() if p.name == "w"]
                        if len(wj_proc): del hists[wj_proc[0]]
                    hists[proc][...] = tmp_fakes

            output[config] = hists
        return output

    def flatten_dy(task, hists, category_inst=None, **kwargs):
        if not hists:
            return hists
        for the_reg, h_single_reg in hists.items():
            if ("fake" in the_reg) or ("gtau" in the_reg):
                continue
            dy_procs = [p for p in h_single_reg.keys() if p.is_mc and "dy" in p.name]
            for p in dy_procs:
                dy_hist = h_single_reg[p].copy()
                if not dy_hist.empty():
                    mean_val = np.average(dy_hist.view().value, axis=1)
                    variance = np.average(dy_hist.view().variance, axis=1) / dy_hist.shape[-1]
                    hists[the_reg][p].view().value = mean_val
                    hists[the_reg][p].view().variance = variance
                else:
                    print(f"DY histogram is empty for {the_reg}")
        return hists

    def symmetrize_signal(task, hists, category_inst=None, **kwargs):
        if not hists:
            return hists
        for the_reg, h_single_reg in hists.items():
            if ("fake" in the_reg) or ("gtau" in the_reg):
                continue
            signal_procs = [p for p in h_single_reg.keys() if p.is_mc and p.has_tag("signal")]
            for p in signal_procs:
                the_hist = h_single_reg[p].copy()
                if the_hist.empty():
                    print(f"signal histogram is empty for {the_reg}")
                    continue

                n_bins = the_hist.shape[-1]

                def get_val(idx, err=None):
                    the_field = "variance" if err else "value"
                    return getattr(hists[the_reg][p].view(), the_field)[:, idx]

                def set_val(hists, idx, val, var):
                    hists[the_reg][p].view().value[:, idx] = val
                    hists[the_reg][p].view().variance[:, idx] = var
                    return hists

                if ("htt_cpo" in p.name) or ("htt_sm" in p.name):
                    for idx in range(n_bins // 2):
                        idxs = [idx, n_bins - idx - 1]
                        val = np.average([get_val(idy) for idy in idxs], axis=0)
                        var = np.average([get_val(idy, err=True) for idy in idxs], axis=0) / 2.0
                        hists = set_val(hists, idxs[0], val, var)
                        hists = set_val(hists, idxs[1], val, var)
                elif "htt_mm" in p.name:
                    for idx in range(n_bins // 4):
                        idxs = [idx, (n_bins // 2) - idx - 1]
                        val = np.average([get_val(idy) for idy in idxs], axis=0)
                        var = np.average([get_val(idy, err=True) for idy in idxs], axis=0) / 2.0
                        hists = set_val(hists, idxs[0], val, var)
                        hists = set_val(hists, idxs[1], val, var)

                        idxs = [idy + (n_bins // 2) for idy in idxs]
                        val = np.average([get_val(idy) for idy in idxs], axis=0)
                        var = np.average([get_val(idy, err=True) for idy in idxs], axis=0) / 2.0
                        hists = set_val(hists, idxs[0], val, var)
                        hists = set_val(hists, idxs[1], val, var)
                else:
                    raise NotImplementedError(
                        "Signal does not have any of these tags [sm, mm, cpo]. "
                        "I don't know how to symmetrize it."
                    )
        return hists

    def blind_sr(task, inputs, variable_name=None, category_name=None, **kwargs):
        outputs = {}
        for config, hists in inputs.items():
            sr_cats = [c for c in config.categories.names() if "_sr" in c]
            out_h = hists.copy()
            shifts = getattr(task, "shifts", ["nominal"])

            for p, h in hists.items():
                if "data" not in p.name:
                    continue
                tmp_arr = h.view()
                for the_cat in sr_cats:
                    for shift_name in shifts:
                        tmp_arr[find_idxs(h, the_cat, shift_name)].value = 1e-2
                        tmp_arr[find_idxs(h, the_cat, shift_name)].variance = 0
                out_h[p][...] = tmp_arr
            outputs[config] = out_h
        return outputs

    def add_qcd_and_wj(task, hists, category_inst=None, **kwargs):
        qcd_proc = _get_config_process(task.config_inst, "qcd")
        if qcd_proc is None:
            logger.warning("task config has no process 'qcd'; skipping add_qcd_and_wj")
            return hists

        fake_hists = []
        for proc, h in hists.items():
            is_fake_proc = ("qcd" in proc.name) or ("wj" in proc.name)
            if is_fake_proc:
                fake_hists.append(h)

        if not fake_hists:
            return hists

        fake_hist = sum(fake_hists[1:], fake_hists[0].copy())
        hists[qcd_proc] = fake_hist
        return hists

    def ensure_zl_hist(task, hists, category_inst=None, **kwargs):
        has_dy_z2ll = [
            the_proc
            for the_proc in hists.keys()
            if ("dy_z2mumu" in the_proc.name) or ("dy_z2ee" in the_proc.name)
        ]
        from cmsdb.processes.ewk import dy_z2ee, dy_z2tautau, dy_z2mumu

        if len(has_dy_z2ll) == 0:
            tmp_hists = [the_h for the_proc, the_h in hists.items() if "dy_z2tautau" in the_proc.name]
            if len(tmp_hists):
                tmp_h = tmp_hists[0]
                hists[dy_z2mumu] = tmp_h.copy().reset()
                hists[dy_z2ee] = tmp_h.copy().reset()
        return hists

    def make_inclusive_hists(task, inputs, variable_name=None, category_name=None, **kwargs):
        outputs = {}
        for config, hists in inputs.items():
            bkg_type = ""

            print(f"Using {bkg_type} subcategories to make inclusive hists.")
            for incl in task.categories:
                sr = config.get_category(incl)
                if 'ff_method' in task.hist_hooks:
                    sr.label = sr.label.replace("prompt lep.", f"QCD Wj and tt data-driven")
                elif 'qcd' in task.hist_hooks:
                    sr.label = sr.label.replace("fake incl.", f"QCD data-driven")
                else:
                    pass
            
            decay_ch = ["tau2pi", "tau2rho", "tau2a1", "tau2a1_3pr"]
            out_h = hists.copy()
            shifts = getattr(task, "shifts", ["nominal"])
            for shift_name in shifts:
                for incl in task.categories:
                    print(f"Preparing hists for {incl} {shift_name} from {config.name}")
                    for p, h in hists.items():
                        tmp_arr = h.view()
                        subhists = []
                        for the_decay in decay_ch:
                            cat_name = "__".join((incl, the_decay))
                            print(f"Adding :{cat_name} for {p.name} from {config.name}")
                            loc_dict = {"category": hist.loc(cat_name), "shift": hist.loc(shift_name)}
                            
                            if cat_name in h.axes['category'] and shift_name in h.axes['shift']:
                                subhists.append(h[loc_dict])
                            else:
                                print(f"WARNING: didn't find {cat_name} for {p.name} from {config.name}")
                        
                        if incl in h.axes['category'] and shift_name in h.axes['shift']:
                            if len(subhists):
                                incl_h = sum(subhists[1:], subhists[0].copy())
                                tmp_arr[find_idxs(h, incl, shift_name)].value = incl_h.view().value
                                tmp_arr[find_idxs(h, incl, shift_name)].variance = incl_h.view().variance
                            else:
                                tmp_arr[find_idxs(h, incl, shift_name)].value = 0
                                tmp_arr[find_idxs(h, incl, shift_name)].variance = 0
                        out_h[p][...] = tmp_arr
            outputs[config] = out_h
        return outputs

    def add_cats(task, inputs, variable_name=None, category_name=None, **kwargs):
        outputs = {}
        for config, hists in inputs.items():
            incl = "cat_mutau_sr"
            subcats = ["cat_mutau_dr_num_wj__prompt", "cat_mutau_dr_num_wj__jet_fakes"]
            out_h = hists.copy()
            shifts = getattr(task, "shifts", "nominal")
            for p, h in hists.items():
                tmp_arr = h.view()
                subhists = []
                for shift in shifts:
                    for subcat in subcats:
                        loc_dict = {"category": hist.loc(subcat), "shift": hist.loc(shift)}
                        subhists.append(h[loc_dict])
                    incl_h = sum(subhists[1:], subhists[0].copy())
                    tmp_arr[find_idxs(h, incl, shift)].value = incl_h.view().value
                    tmp_arr[find_idxs(h, incl, shift)].variance = incl_h.view().variance
                out_h[p][...] = tmp_arr
            outputs[config] = out_h
        return outputs

    def order_hists(task, inputs, variable_name=None, category_name=None, **kwargs):
        outputs = {}
        for config, hists in inputs.items():
            cat_name = task.categories
            cat = config.get_category(cat_name[0])
            #cat.label = "Signal Region"
            out_hists = {}
            procs = [
                p
                for p in hists.keys()
                if (p.name != "qcd") and (p.name != "dy_tt_m50") and (p.name != "dy_ll_m50")
            ]
            qcd = [p for p in hists.keys() if p.name == "qcd"]
            dy_tt = [p for p in hists.keys() if p.name == "dy_tt_m50"]
            dy_ll = [p for p in hists.keys() if p.name == "dy_ll_m50"]

            if len(qcd):
                out_hists[qcd[0]] = hists[qcd[0]]
            for p in procs:
                out_hists[p] = hists[p]
            if len(dy_ll):
                out_hists[dy_ll[0]] = hists[dy_ll[0]]
            if len(dy_tt):
                out_hists[dy_tt[0]] = hists[dy_tt[0]]

            outputs[config] = out_hists
        return outputs

    analysis.x.hist_hooks = {
        "qcd": qcd_estimation,
        "add_qcd_and_wj": add_qcd_and_wj,
        "ff_method": ff_method,
        "closure_test": ff_closure_test,
        "flatten_dy": flatten_dy,
        "symmetrize_signal": symmetrize_signal,
        "blind_sr": blind_sr,
        "ensure_zl_hist": ensure_zl_hist,
        "order": order_hists,
        "incl": make_inclusive_hists,
        "add_cats": add_cats,
    }