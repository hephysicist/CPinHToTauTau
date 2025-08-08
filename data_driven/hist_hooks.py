"""
Histogram hooks.
"""

from collections import defaultdict

import law
import order as od

import scinum as sn

from columnflow.util import maybe_import
from law.util import DotDict

np = maybe_import("numpy")
ak = maybe_import("awkward")
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

def get_signal_hists(hists: dict)-> hist.Hist :
    hists = [h 
             for p, h in hists.items() 
             if p.is_mc and p.has_tag("signal")]
    return hists

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
    flat_tf = True
    
    
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
        
        
        from cmsdb.processes.qcd import qcd
        hists_sr = hists[sr.name].copy()
        h_donor_name = list(hists_sr.keys())[0]
        hists_sr[qcd] = hists_sr[h_donor_name].copy().reset()
        
        if not data_ar.empty():
            if flat_tf:
                num = ak.sum(data_num.values() - mc_num.values())
                den = ak.sum(data_den.values() - mc_den.values())

                if (num > 0) and (den > 0):
                    tf = num/den
                else:
                    tf = 1. 
            
            else:
                data_num.values() - mc_num.values()
                data_den.values() - mc_den.values()
            
        
                mask = ((num > 0) & (den > 0))
            
                tf = num/den
                tf = ak.where((num>0) & (den>0), tf, np.ones_like(num))
        
                tf_err2 = ((np.sum(data_num.variances()) + np.sum(mc_num.variances()))/den**2 + 
                    tf**2/den**2 *(np.sum(data_den.variances()) + np.sum(mc_den.variances())))
            hists_sr[qcd].view().value = np.maximum(data_ar.values() -  mc_ar.values(), 0.) * tf
            
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
                        if ('cpo' in p.name) or ('sm' in p.name):
                            for idx in range(n_bins//2):
                                idxs = [idx,n_bins-idx-1]
                                val = np.average([get_val(idy) for idy in idxs])
                                var = np.average([get_val(idy,err=True) for idy in idxs])/2.
                                hists = set_val(hists, idxs[0], val, var)
                                hists = set_val(hists, idxs[1], val, var)
                        elif ('mm' in p.name):
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
    
    
    def blind_sr(task, hists, category_inst):
        if not hists:
            return hists
        data_proc = [the_proc for the_proc in hists.keys() if 'data' in the_proc.name]
        hists[data_proc[0]].view().value[:,1:] = 0
        hists[data_proc[0]].view().variance[:,1:] = 0
        return hists
    
    def add_qcd_and_wj(task, hists, category_inst):
        fake_hists = []
        for (proc, h) in hists.items():
            is_fake_proc = ('qcd' in proc.name) or ('wj' in proc.name )
            if is_fake_proc: fake_hists.append(h)
        fake_hist = sum(fake_hists[1:], fake_hists[0].copy())
        from cmsdb.processes.qcd import qcd
        hists[qcd] = fake_hist
        return hists
    
    
    config.x.hist_hooks = {
        "good_old_abcd"             : qcd_estimation,
        "add_qcd_and_wj"            : add_qcd_and_wj,
        "ff_method"                 : ff_method,
        "ff_method_dr_closure_test" : ff_method_dr_closure_test,
        "flatten_dy"                : flatten_dy,
        "symmetrize_signal"         : symmetrize_signal,
        "blind_sr"                  : blind_sr,
    }