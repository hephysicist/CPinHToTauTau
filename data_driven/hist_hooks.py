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




def calc_yields(hists: dict)-> hist.Hist :
    data_hists = [h for p, h in hists.items() if p.is_data]
    data_hist = sum(data_hists[1:], data_hists[0].copy())
    data = data_hist.values()
    
    mc_hists = [h for p, h in hists.items() if (p.is_mc and not p.has_tag("signal"))]
    mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
    mc = mc_hist.values()
    
    wj_hists = [h for p, h in hists.items() if (p.is_mc and (not p.has_tag("signal")) and 'wj' in p.name)]
    wj_hist = sum(wj_hists[1:], wj_hists[0].copy())
    wj = wj_hist.values()
    
    wj_ratio =  wj / np.maximum(data, 1)
    wj_err   =  np.where(data > 0,
                            rel_err(h_arr=[wj_hist,data_hist]) * wj_ratio,
                            np.ones_like(wj_ratio)) 
    
    qcd_ratio =  np.maximum((data - mc), 0)/ np.maximum(data, 1)
    qcd_err   =  np.where((data - mc) > 0,
                        rel_err(h_arr=[wj_hist,data_hist]) * qcd_ratio,
                        np.ones_like(wj_ratio)) 
    
    return {'wj': wj_ratio[0],
            'wj_err': wj_err[0],
            'qcd': qcd_ratio[0],
            'qcd_err': qcd_err[0],}

def get_data_hist(hists: dict)-> hist.Hist :
    data_hists = [h for p, h in hists.items() if p.is_data]
    data_hist = sum(data_hists[1:], data_hists[0].copy())
    return data_hist

def get_mc_hist(hists: dict, remove_wj=False)-> hist.Hist :
    data_hists = [h for 
                  p, h in hists.items() 
                  if p.is_mc 
                  and not p.has_tag("signal")
                  and ((remove_wj * ('wj' not in p.name)) or ~remove_wj) ]
    data_hist = sum(data_hists[1:], data_hists[0].copy())
    return data_hist

def rel_err(h_arr=[], err_arr=[]):
    if len(h_arr):  sum_var = np.zeros_like(h_arr[0].values())
    else: sum_var = np.zeros_like(err_arr[0])
    for x in h_arr: sum_var += x.variances()/np.maximum(x.values()**2, 1)
    for the_arr in err_arr: sum_var += err_arr
    return np.sqrt(sum_var)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    
    
    def qcd_estimation(task, hists, category_inst):
        
        def get_hists_from_reg(config: od.Config, hists: dict, region: str)-> hist.Hist :
            hists_ = hists[region]
            cat_id = config.get_category(region).id
            data_hists = []
            mc_hists = []
            for (proc, h) in hists_.items():
                if proc.is_data:
                    data_hists.append(h)
                elif proc.is_mc and not proc.has_tag("signal"):
                    mc_hists.append(h)
            
            mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
            data_hist = sum(data_hists[1:], data_hists[0].copy())
            
            return data_hist, mc_hist
        
        sr = category_inst
        data_num, mc_num = get_hists_from_reg(config, hists,sr.aux['abcd_regs']['dr_num'])
        data_den, mc_den = get_hists_from_reg(config, hists, sr.aux['abcd_regs']['dr_den']) 
        data_ar, mc_ar = get_hists_from_reg(config, hists,sr.aux['abcd_regs']['ar'])
        num = np.sum(data_num.values()) - np.sum(mc_num.values())
        den = np.sum(data_den.values()) - np.sum(mc_den.values()) 
        
        tf = num/np.maximum(den,0.001)
        tf_err2 = ((np.sum(data_num.variances()) + np.sum(mc_num.variances()))/den**2 + 
                  tf**2/den**2 *(np.sum(data_den.variances()) + np.sum(mc_den.variances())))
        
        from cmsdb.processes.qcd import qcd
        hists_sr = hists[sr.name].copy()
        h_donor_name = list(hists_sr.keys())[0]
        hists_sr[qcd] = hists_sr[h_donor_name].copy().reset()
        hists_sr[qcd].view().value = np.maximum(data_ar.values() -  mc_ar.values(), 0.) * tf
        #hists_sr[qcd].view().variance = data_ar.values()**2 * tf_err2 + data_ar.variances() * (tf**2)
        return hists_sr
    
    def ff_method(task, hists, category_inst):
        if not hists:
            return hists
        
        sr = category_inst
        
        hist_qcd = get_data_hist(hists[sr.aux['ff_regs']['ar_qcd']])
        hist_wj  = get_data_hist(hists[sr.aux['ff_regs']['ar_wj']])
        yields = calc_yields(hists[sr.aux['ff_regs']['ar_yields']])
      
        fakes_wj  = hist_wj.values() * yields['wj']
        fakes_qcd = hist_qcd.values() * yields['qcd']
        
        hists_sr = hists[sr.name].copy()
        
        #Remove wj histogram from the signal region set
        from cmsdb.processes.qcd import jet_fakes,qcd
        
        wj_proc = [p for p in hists_sr.keys() if 'wj' in p.name]
        
        tmp_h = list(hists_sr.values())
        tmp_h = sum(tmp_h[1:],tmp_h[0].copy())
        hists_sr[jet_fakes] = tmp_h.copy().reset()
        hists_sr[jet_fakes].view().value = fakes_wj
        
        hists_sr[qcd] = tmp_h.copy().reset()
        hists_sr[qcd].view().value = fakes_qcd
        
        del hists_sr[wj_proc[0]]   
        return hists_sr
    
    def ff_method_dr_closure_test(task, hists, category_inst):
        if not hists:
            return hists
        cat_tag = category_inst.name.split('__')
        if len(cat_tag) > 1:
            cat_tag = '__' + cat_tag[1] 
        else:
            cat_tag = ''
        sr = task.config_inst.get_category('cat_mutau_sr' + cat_tag)
        
        hist_qcd = get_data_hist(hists[sr.aux['ff_regs']['dr_den_qcd_w_ff']])
        hist_wj  = get_data_hist(hists[sr.aux['ff_regs']['dr_den_wj_w_ff']])
        
        hist_qcd_mc = get_mc_hist(hists[sr.aux['ff_regs']['dr_den_qcd_w_ff']])
        hist_wj_mc  = get_mc_hist(hists[sr.aux['ff_regs']['dr_den_wj_w_ff']])

        fakes_wj  = (hist_wj.values()  - hist_wj_mc.values()) 
        fakes_qcd = (hist_qcd.values() - hist_qcd_mc.values()) 
        
        hists_sr = hists[category_inst.name].copy()
        tmp_h = list(hists_sr.values())
        tmp_h = sum(tmp_h[1:],tmp_h[0].copy())
        #Remove wj histogram from the signal region set
        from cmsdb.processes.qcd import jet_fakes,qcd
        wj_proc = [p for p in hists_sr.keys() if 'wj' in p.name]
        if 'wj' in category_inst.name:
            hists_sr[jet_fakes] = tmp_h.copy().reset()
            hists_sr[jet_fakes].view().value = fakes_wj
            del hists_sr[wj_proc[0]]
        else:
            hists_sr[qcd] = tmp_h.copy().reset()
            hists_sr[qcd].view().value = fakes_qcd
            
        return hists_sr
    

    config.x.hist_hooks = {
        "good_old_abcd"  : qcd_estimation,
        "ff_method" : ff_method,
        "ff_method_dr_closure_test": ff_method_dr_closure_test,
    }
