"""
Histogram hooks.
"""

from collections import defaultdict

import law
import order as od

import scinum as sn

from columnflow.util import maybe_import
from law.util import InsertableDict
import warnings

np = maybe_import("numpy")
ak = maybe_import("awkward")
hist = maybe_import("hist")

logger = law.logger.get_logger(__name__)

def calc_yields(hists: dict, locator: dict, fake_locator: dict)-> hist.Hist :
    #Get data hists
    data_hists = [h for p, h in hists.items() if p.is_data]
    data_hist = sum(data_hists[1:], data_hists[0].copy())
    data = data_hist[locator].values()
    
    #Get all monte-carlo hists
    mc_hists = [h for p, h in hists.items() if (p.is_mc and not p.has_tag("signal"))]
    mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
    mc = mc_hist[locator].values()
    mc_fakes = mc_hist[fake_locator].values()

    wj_ratio =  (mc_fakes) / np.maximum(data - mc , 1)
    wj_err   = np.ones_like(wj_ratio)
    
    qcd_ratio =  np.maximum((data - (mc_fakes + mc)), 0)/ np.maximum(data - mc, 1)
    qcd_err   = np.ones_like(wj_ratio) 
    # np.where(data > 0,
    #                         rel_err(h_arr=[wj_hist[fake_locator],data_hist[locator]]) * wj_ratio,
    #                         np.ones_like(wj_ratio)) 
    
   
    # np.where((data - mc) > 0,
    #                     rel_err(h_arr=[wj_hist[fake_locator],data_hist[locator]]) * qcd_ratio,
    #                     np.ones_like(wj_ratio)) 
    
    return {'wj': wj_ratio,
            'wj_err': wj_err,
            'qcd': qcd_ratio,
            'qcd_err': qcd_err,}

def get_data_hist(hists: dict)-> hist.Hist :
    data_hists = [h for p, h in hists.items() if p.is_data]
    if len(data_hists): return sum(data_hists[1:], data_hists[0].copy())
    else: return None

def get_mc_hist(hists: dict, remove_wj=False)-> hist.Hist :
    data_hists = [h for 
                  p, h in hists.items() 
                  if p.is_mc 
                  and not p.has_tag("signal")
                  and ((remove_wj * ('wj' not in p.name)) or not remove_wj) ]
    if len(data_hists): return sum(data_hists[1:], data_hists[0].copy())
    else: return None
    

def get_signal_hists(hists: dict)-> hist.Hist :
    hists = [h 
             for p, h in hists.items() 
             if p.is_mc and p.has_tag("signal") and (p.name != 'qcd')]
    return hists

def rel_err(h_arr=[], err_arr=[]):
    if len(h_arr):  sum_var = np.zeros_like(h_arr[0].values())
    else: sum_var = np.zeros_like(err_arr[0])
    for x in h_arr: sum_var += x.variances()/np.maximum(x.values()**2, 1)
    for the_arr in err_arr: sum_var += err_arr
    return np.sqrt(sum_var)

def find_idxs(h: hist.Hist, cat: str, shift: str):
    return h.axes[:-1].index(cat,shift)

def add_hist_hooks(analysis: od.Analysis) -> None:
    """
    Add histogram hooks to a configuration.
    """
    flat_tf = True
    def qcd_estimation(task, inputs): #cf0p3
        output = {}
        for config, hists in inputs.items():
            sr_cats = [c for c in config.categories.names() if '_sr' in c]
            for the_cat in sr_cats:
                print(f'producing qcd for {the_cat}')
                sr = config.get_category(the_cat)
                d = {}
                mc = {}
                cr_cat = ''
                for reg_name, full_name in sr.aux['abcd_regs'].items():   
                    loc_dict = {'category': hist.loc(full_name),'shift': hist.loc(task.shift)}
                    full_d = get_data_hist(hists)
                    full_mc = get_mc_hist(hists)   
                    if full_name in list(full_d.axes[0]): d[reg_name] = get_data_hist(hists)[loc_dict]
                    if full_name in list(full_mc.axes[0]): mc[reg_name] = get_mc_hist(hists)[loc_dict]
                    if reg_name == 'dr_num': cr_cat = full_name
                        
                from cmsdb.processes.qcd import qcd
                h_donor_name  = list(hists.keys())[0]
                if qcd not in hists.keys():
                    hists[qcd] = hists[h_donor_name].copy().reset()
                    tmp_arr = hists[qcd].view()
                if 'ar' in d.keys() and not d['ar'].empty():
                    if flat_tf:
                        tf = 1
                        #num = ak.sum(d['dr_num'].values() - mc['dr_num'].values())
                        #den = ak.sum(d['dr_den'].values() - mc['dr_den'].values())
                        #if (num > 0) and (den > 0):
                        #    tf = num/den
                        #else:
                        #    tf = 1. 
                    
                    else:
                        num = d['dr_num'].values() - mc['dr_num'].values()
                        den = d['dr_den'].values() - mc['dr_den'].values()

                        # mask = ((num > 0) & (den > 0))
                    
                        # tf = num/den
                        # tf = ak.where((num>0) & (den>0), tf, np.ones_like(num))
                
                        # tf_err2 = ((np.sum(d['dr_num'].variances()) + np.sum(mc['dr_num'].variances()))/den**2 + 
                        #     tf**2/den**2 *(np.sum(d['dr_den'].variances()) + np.sum(mc['dr_den'].variances())))
                    
                    val = np.maximum(d['ar'].values() - mc['ar'].values(), 0.) * tf
                    var = (d['ar'].view().variance + mc['ar'].view().variance) * tf**2
                  
                    tmp_arr[find_idxs(hists[qcd], the_cat, task.shift)].value = val
                    tmp_arr[find_idxs(hists[qcd], the_cat, task.shift)].variance = var
                    if len(cr_cat):
                        if ('dr_den' in d.keys()) and ('dr_den' in mc.keys()):  
                            cr_val = np.maximum(d['dr_den'].values() - mc['dr_den'].values(), 0)
                            cr_var = d['dr_den'].variances() - mc['dr_den'].variances()
                        elif ('dr_den' in d.keys()):
                            cr_val = d['dr_den'].values() 
                            cr_var = d['dr_den'].variances()
                        else:
                            cr_val = 0
                            cr_var = 0
                        tmp_arr[find_idxs(hists[qcd], cr_cat, task.shift)].value = cr_val
                        tmp_arr[find_idxs(hists[qcd], cr_cat, task.shift)].variance = cr_var
                else:
                    print("*** WARNING: AR data histogam doesn't exist or empty! ***")
                    
            hists[qcd][...] = tmp_arr
            output[config] = hists
            return output
    
    def ff_method(task, inputs): #cf0p3
        from cmsdb.processes.qcd import jet_fakes,qcd
        output = {}
        for config, hists in inputs.items():
            sr_cats = [c for c in config.categories.names() if ('_sr' in c) and ('no_mt' not in c) ]
            for the_cat in sr_cats:
                print(f'Applying fake factor method on {the_cat}')
                sr = config.get_category(the_cat)
                data = get_data_hist(hists)
                mc = get_mc_hist(hists)
               
                locator = lambda reg:  {'category': hist.loc(sr.aux['ff_regs'][reg]),'shift': hist.loc(task.shift)}
                
                data_qcd = data[locator('ar_qcd')]
                mc_qcd = mc[locator('ar_qcd')]
                data_wj = data[locator('ar_wj')]
                mc_wj = mc[locator('ar_wj')]
                #from IPython import embed; embed()
                yields = calc_yields(hists, locator('ar_yields'), locator('ar_yields_fakes'))
                
                
                h_donor_name  = list(hists.keys())[0]
                if qcd not in hists.keys():
                    hists[qcd] = hists[h_donor_name].copy().reset()
                    tmp_qcd = hists[qcd].view()
                
                if jet_fakes not in hists.keys():
                    hists[jet_fakes] = hists[h_donor_name].copy().reset()
                    tmp_fakes = hists[jet_fakes].view()
                    
                tmp_qcd[find_idxs(hists[qcd], the_cat, task.shift)].value = np.maximum(data_qcd.values() - mc_qcd.values(), 0) * yields['qcd']
                tmp_qcd[find_idxs(hists[qcd], the_cat, task.shift)].variance = data_qcd.variances() + mc_qcd.variances()
                
                tmp_fakes[find_idxs(hists[jet_fakes], the_cat, task.shift)].value = np.maximum(data_wj.values() - mc_wj.values(), 0) * yields['wj']
                tmp_fakes[find_idxs(hists[jet_fakes], the_cat, task.shift)].variance = data_wj.variances() + mc_wj.variances()
                
            hists[qcd][...] = tmp_qcd
            hists[jet_fakes][...] = tmp_fakes
            
            [wj_proc] = [p for p in hists.keys() if 'wj' in p.name]
            del hists[wj_proc]   
            output[config] = hists
            return output
    
    def ff_closure_test(task, inputs): #cf0p3
        from cmsdb.processes.qcd import jet_fakes,qcd
        output = {}
        for config, hists in inputs.items():
            dr_num_cats = task.categories
            for the_cat in dr_num_cats:
                print(f'Performin closure test for {the_cat}')
                sr_name = the_cat.replace('dr_num_wj','sr').replace('dr_num_qcd','sr')
                sr = config.get_category(sr_name)
                data = get_data_hist(hists)
                mc = get_mc_hist(hists)
                locator = lambda reg:  {'category': hist.loc(sr.aux['ff_regs'][reg]),'shift': hist.loc(task.shift)}
                iter_dict = {'wj': jet_fakes, 'qcd':qcd}
                h_donor_name  = list(hists.keys())[0]
                for name, proc in iter_dict.items():
                    if proc not in hists.keys():
                        hists[proc] = hists[h_donor_name].copy().reset()
                        tmp_fakes = hists[proc].view()
                        h_d = data[locator(f'dr_den_{name}_w_ff')]
                        h_mc = mc[locator(f'dr_den_{name}_w_ff')]
                        num = sr.aux['ff_regs'][f'dr_num_{name}']
                        tmp_fakes[find_idxs(hists[jet_fakes], num, task.shift)].value = np.maximum(h_d.values() - h_mc.values(), 0)
                        tmp_fakes[find_idxs(hists[jet_fakes], num, task.shift)].variance = h_d.variances() + h_mc.variances()
                        hists[proc][...] = tmp_fakes
                    if name == 'wj':
                        p_list = [p for p in hists.keys() if 'wj' in p.name]
                        if len(p_list):
                            del hists[p_list[0]]          
            output[config] = hists
            return output
    
    # def ff_method_dr_closure_test(task, hists):
    #     if not hists:
    #         return hists
    #     cat_tag = category_inst.name.split('__')
    #     if len(cat_tag) > 1:
    #         cat_tag = '__' + cat_tag[1] 
    #     else:
    #         cat_tag = ''
    #     sr = task.config_inst.get_category('cat_mutau_sr' + cat_tag)
        
    #     hist_qcd = get_data_hist(hists[sr.aux['ff_regs']['dr_den_qcd_w_ff']])
    #     hist_wj  = get_data_hist(hists[sr.aux['ff_regs']['dr_den_wj_w_ff']])
        
    #     hist_qcd_mc = get_mc_hist(hists[sr.aux['ff_regs']['dr_den_qcd_w_ff']])
    #     hist_wj_mc  = get_mc_hist(hists[sr.aux['ff_regs']['dr_den_wj_w_ff']])

    #     fakes_wj  = (hist_wj.values()  - hist_wj_mc.values()) 
    #     fakes_qcd = (hist_qcd.values() - hist_qcd_mc.values()) 
        
    #     hists_sr = hists[category_inst.name].copy()
    #     tmp_h = list(hists_sr.values())
    #     tmp_h = sum(tmp_h[1:],tmp_h[0].copy())
    #     #Remove wj histogram from the signal region set
    #     from cmsdb.processes.qcd import jet_fakes,qcd
    #     wj_proc = [p for p in hists_sr.keys() if 'wj' in p.name]
    #     if 'wj' in category_inst.name:
    #         hists_sr[jet_fakes] = tmp_h.copy().reset()
    #         hists_sr[jet_fakes].view().value = fakes_wj
    #         del hists_sr[wj_proc[0]]
    #     else:
    #         hists_sr[qcd] = tmp_h.copy().reset()
    #         hists_sr[qcd].view().value = fakes_qcd
            
    #     return hists_sr
    
    def flatten_dy(task, hists, category_inst):
        if not hists:
            return hists
        for the_reg, h_single_reg in hists.items():
            if ('fake' in the_reg) or ('gtau' in the_reg):
                pass
            else:
                dy_procs = [p for p in h_single_reg.keys() if p.is_mc and 'dy' in p.name]
                for p in dy_procs:
                    dy_hist = h_single_reg[p].copy()
                    if not dy_hist.empty():
                        mean_val = np.average(dy_hist.view().value, axis=1)
                        variance =np.average(dy_hist.view().variance, axis=1)/dy_hist.shape[-1]
                        #print(f"Perform dy flattening for {the_reg}")
                        #print(f"Before: {hists[the_reg][p].view().value}")
                        hists[the_reg][p].view().value = mean_val
                        hists[the_reg][p].view().variance = variance
                        #print(f"After: {hists[the_reg][p].view().value}")
                    else:
                        print(f"DY histogram is empty for {the_reg}")
        return hists

    def symmetrize_signal(task, hists, category_inst):
        if not hists:
            return hists
        for the_reg, h_single_reg in hists.items():
            if ('fake' in the_reg) or ('gtau' in the_reg):
                pass
            else:
                signal_procs = [p for p in h_single_reg.keys() if p.is_mc and p.has_tag("signal")]
                for p in signal_procs:
                    the_hist = h_single_reg[p].copy()
                    if not the_hist.empty():
                        n_bins = the_hist.shape[-1]
                        def get_val(idx, err=None):
                            the_field = 'variance' if err else 'value'
                            return getattr(hists[the_reg][p].view(), the_field)[:,idx]
                        def set_val(hists, idx, val, var):
                            hists[the_reg][p].view().value[:,idx]= val
                            hists[the_reg][p].view().variance[:,idx]= var
                            return hists 
                        if ('htt_cpo' in p.name) or ('htt_sm' in p.name):
                            for idx in range(n_bins//2):
                                idxs = [idx,n_bins-idx-1]
                                val = np.average([get_val(idy) for idy in idxs])
                                var = np.average([get_val(idy,err=True) for idy in idxs])/2.
                                hists = set_val(hists, idxs[0], val, var)
                                hists = set_val(hists, idxs[1], val, var)
                        elif ('htt_mm' in p.name):
                            for idx in range(n_bins//4):
                                #left part of the distribution 
                                idxs = [idx,(n_bins//2)-idx-1]
                                val = np.average([get_val(idy) for idy in idxs])
                                var = np.average([get_val(idy,err=True) for idy in idxs])/2.
                                hists = set_val(hists, idxs[0], val, var)
                                hists = set_val(hists, idxs[1], val, var)
                                #right part of the distribution
                                idxs = [idy+(n_bins//2) for idy in idxs]
                                val = np.average([get_val(idy) for idy in idxs])
                                var = np.average([get_val(idy,err=True) for idy in idxs])/2.
                                hists = set_val(hists, idxs[0], val, var)
                                hists = set_val(hists, idxs[1], val, var)
                        else:
                            raise NotImplementedError(f"Signal does not have any of this tags [sm, mm, cpo]. I don't know how to symmetrize it.")
                        norm_after = np.sum(hists[the_reg][p].view().value,axis=1)
                    else:
                        print(f"signal histogram is empty for {the_reg}")
        return hists


    def blind_sr(task, inputs):
        ouputs = {}
       
        for config, hists in inputs.items():
            sr_cats = [c for c in config.categories.names() if '_sr' in c]
            out_h = hists.copy()
            for p, h in hists.items():
                if 'data' in p.name:
                    tmp_arr = h.view()
                    for the_cat in sr_cats:
                        loc_dict = {'category': hist.loc(the_cat),'shift': hist.loc(task.shift)}
                        tmp_arr[find_idxs(h, the_cat, task.shift)].value[-5:] = 0
                        tmp_arr[find_idxs(h, the_cat, task.shift)].variance[-5:] = 0

                    out_h[p][...] = tmp_arr
            ouputs[config] = out_h
        return ouputs

    def add_qcd_and_wj(task, hists, category_inst):
        fake_hists = []
        for (proc, h) in hists.items():
            is_fake_proc = ('qcd' in proc.name) or ('wj' in proc.name )
            if is_fake_proc: fake_hists.append(h)
        fake_hist = sum(fake_hists[1:], fake_hists[0].copy())
        from cmsdb.processes.qcd import qcd
        hists[qcd] = fake_hist
        return hists

    def ensure_zl_hist(task, hists, category_inst):
        has_dy_z2ll = [the_proc 
                       for the_proc in hists.keys() 
                       if ('dy_z2mumu' in the_proc.name) or ('dy_z2ee' in the_proc.name)]
        from cmsdb.processes.ewk import dy_z2ee,dy_z2tautau
        if len(has_dy_z2ll) == 0:
            tmp_h = [the_h for proc , the_h in hists if ('dy_z2tautau' in the_proc.name)][0] 
            hists[dy_z2mumu] = tmp_h.copy().reset()
            hists[dy_z2ee] = tmp_h.copy().reset()
        return hists
    
    
    def make_inclusive_hists(task, inputs): #cf0p3
        
        ouputs = {}
        for config, hists in inputs.items():
            incl = 'cat_mutau_sr'
            bkg_type = ['prompt','jet_fakes']
            hig_cats = ['hig__cat0','hig__cat1', 'hig__cat2'] 
            decay_ch = ['tau2pi','tau2rho','tau2a1','tau2a1_3pr']
            out_h = hists.copy()
            for p, h in hists.items():
                tmp_arr = h.view()
                subhists = []
                for bkg in bkg_type:
                    if p.is_data and bkg=='jet_fakes': 
                        print('Not adding data second time')
                        continue #Otherwise we double count data from prompt cats and from the jet_fakes cats
                    for the_cat in hig_cats:
                        for the_decay in decay_ch:
                            cat_name = '__'.join((incl,bkg,the_cat,the_decay))
                            print(f'Adding :{cat_name} for {p.name} from {config.name}')
                            loc_dict = {'category': hist.loc(cat_name),
                                        'shift': hist.loc(task.shift)}
                            if cat_name in h.axes[0]:
                                subhists.append(h[loc_dict])
                            else:
                                print(f"WARNING: didn't find {cat_name} for {p.name} from {config.name}")
                incl_h = sum(subhists)
                tmp_arr[find_idxs(h, incl, task.shift)].value = incl_h.view().value
                tmp_arr[find_idxs(h, incl, task.shift)].variance = incl_h.view().variance
                out_h[p][...] = tmp_arr
            ouputs[config] = out_h
        return ouputs

    def add_cats(task, inputs): #cf0p3
        ouputs = {}
        for config, hists in inputs.items():
            incl = 'cat_mutau_sr'
            #subcats = ['cat_mutau_sr','cat_mutau_abcd_dr_num'] #to plot iso inclusive
            subcats = ['cat_mutau_sr','cat_mutau_dr_num_wj'] #to plot mt inclusive
            out_h = hists.copy()
            for p, h in hists.items():
                tmp_arr = h.view()
                subhists = []
                for subcat in subcats:
                    loc_dict = {'category': hist.loc(subcat),'shift': hist.loc(task.shift)}
                    subhists.append(h[loc_dict])
                incl_h = sum(subhists)
                tmp_arr[find_idxs(h, incl, task.shift)].value = incl_h.view().value
                tmp_arr[find_idxs(h, incl, task.shift)].variance = incl_h.view().variance
                out_h[p][...] = tmp_arr
            ouputs[config] = out_h
        return ouputs
    

    def order_hists(task, inputs): #cf0p3
        ouputs = {}
        for config, hists in inputs.items():
            out_hists = {}
            procs = [p for p in hists.keys() if (p.name != 'qcd') and (p.name != 'dy_tt_m50') and (p.name != 'dy_ll_m50')]
            qcd = [p for p in hists.keys() if (p.name == 'qcd')]
            dy_tt = [p for p in hists.keys() if (p.name == 'dy_tt_m50')]
            dy_ll = [p for p in hists.keys() if (p.name == 'dy_ll_m50')]
            if len(qcd): out_hists[qcd[0]] = hists[qcd[0]]
            for p in procs: out_hists[p] = hists[p]
            if len(dy_ll): out_hists[dy_ll[0]] = hists[dy_ll[0]]
            if len(dy_tt): out_hists[dy_tt[0]] = hists[dy_tt[0]]
            ouputs[config] = out_hists
        return ouputs
    
    

    analysis.x.hist_hooks = {
        "qcd"                       : qcd_estimation,
        "add_qcd_and_wj"            : add_qcd_and_wj,
        "ff_method"                 : ff_method,
        "closure_test"              : ff_closure_test,
        "flatten_dy"                : flatten_dy,
        "symmetrize_signal"         : symmetrize_signal,
        "blind_sr"                  : blind_sr,
        "ensure_zl_hist"            : ensure_zl_hist,
        "order"                     : order_hists,
       "incl"                       : make_inclusive_hists,
       "add_cats"                   : add_cats,
       }
