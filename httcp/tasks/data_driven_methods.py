
"""
Task to produce and merge histograms.
"""

import law

from columnflow.tasks.framework.base import Requirements, BaseTask, AnalysisTask, DatasetTask, wrapper_factory
from columnflow.tasks.framework.mixins import (ProducersMixin, 
                                               ChunkedIOMixin, DatasetsProcessesMixin, CategoriesMixin,VariablesMixin
)
from columnflow.tasks.framework.mixins import (
    CalibratorClassesMixin, SelectorClassMixin, ReducerClassMixin, ProducerClassesMixin, HistProducerClassMixin,
    CategoriesMixin, ShiftSourcesMixin, HistHookMixin, MLModelsMixin,
)
from columnflow.tasks.framework.plotting import ProcessPlotSettingMixin

from columnflow.tasks.framework.remote import RemoteWorkflow
from columnflow.tasks.framework.parameters import last_edge_inclusive_inst
from columnflow.tasks.reduction import ReducedEventsUser
from columnflow.tasks.production import ProduceColumns
from columnflow.tasks.histograms import MergeHistograms
from columnflow.util import dev_sandbox, DotDict


class ComputeFakeFactors(
    CalibratorClassesMixin,
    SelectorClassMixin,
    ReducerClassMixin,
    ProducerClassesMixin,
    MLModelsMixin,
    HistProducerClassMixin,
    CategoriesMixin,
    law.LocalWorkflow,
    RemoteWorkflow,
    VariablesMixin,
    DatasetsProcessesMixin,
    law.tasks.RunOnceTask
):
    single_config = False

    sandbox = dev_sandbox(law.config.get("analysis", "default_columnar_sandbox"))

    exclude_index = False

    def store_parts(self) -> law.util.InsertableDict:
        parts = super().store_parts()
        parts.insert_before("version", "datasets", f"datasets_{self.datasets_repr}")
        return parts

    def create_branch_map(self):
        return [
            DotDict({"category": cat_name, "variable": var_name})
            for cat_name in sorted(self.categories)
            for var_name in sorted(self.variables)
        ]

    resolution_task_cls = MergeHistograms

    reqs = Requirements(
        RemoteWorkflow.reqs,
        MergeHistograms=MergeHistograms
    )

    def workflow_requires(self):
        reqs = super().workflow_requires()
        reqs["merged_hists"] = self.requires_from_branch()
        return reqs

    def requires(self):
        req = {}
        for i, config_inst in enumerate(self.config_insts):
            sub_datasets = self.datasets[i]
            req[config_inst.name] = {}
            for d in sub_datasets:
                if d in config_inst.datasets.names():
                    req[config_inst.name][d] = self.reqs.MergeHistograms.req(
                        self,
                        config=config_inst.name,
                        shift='nominal',
                        dataset=d,
                        branch=-1,
                        _exclude={"branches"},
                        _prefer_cli={"variables"},
                    )
        return req
    
    @classmethod
    def req_params(cls, inst: AnalysisTask, **kwargs) -> dict:
        _prefer_cli = law.util.make_set(kwargs.get("_prefer_cli", [])) | {"variables"}
        kwargs["_prefer_cli"] = _prefer_cli
        return super().req_params(inst, **kwargs)
    
    def output(self):
        cfg = self.config_insts[0]
        if len(self.configs) > 1:
            tag = ''
        else:
            tag = cfg.campaign.aux['tag']
        year = cfg.campaign.aux['year']
        channel = cfg.channels.get_first().name
        return {"ff_json": self.target('_'.join(('fake_factors',
                                                 channel,
                                                 str(year),
                                                 tag)) + '.json'),
                "plots": {'_'.join((ff_type,
                                    syst,
                                    f'n_jets_{str(nj)}')): self.target(f"fake_factor_{ff_type}_{syst}_njets_{str(nj)}.png")
                          for syst in ['nominal', 'up', 'down']
                          for ff_type in ['qcd','wj']
                          for nj in [0,1,2]},
                "plots1d": {'_'.join((ff_type,
                                      str(dm),
                                      str(nj))): self.target(f"fake_factor_{ff_type}_PNet_dm_{str(dm)}_njets_{str(nj)}.png")
                          for ff_type in ['qcd','wj']
                          for dm in [0,1,2,10,11]
                          for nj in [0,1,2]},
                "fitres": self.target('_'.join(('fitres',
                                                 channel,
                                                 str(year),
                                                 tag)) + '.json'),
                }

    @law.tasks.RunOnceTask.complete_on_success
    def run(self):
        import hist
        import numpy as np
        from scipy.optimize import curve_fit
        from scipy.special import erf
        import matplotlib.pyplot as plt
        import correctionlib.schemav2 as cs
        from numpy import exp
        plt.figure(dpi=200)
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "monospace",
            "font.monospace": 'Computer Modern Typewriter'
        })
        # preare inputs and outputs
        inputs = self.input()
        outputs = self.output()
        data_hists = None
        mc_hists = None
        
        for the_cfg in self.config_insts:
            cfg_name = the_cfg.name
            print(f'*** config {cfg_name} ***')
            for (dataset_name, the_input) in inputs[cfg_name].items():
                # load input histograms per dataset
                the_hist = the_input["collection"][0]['hists'].targets[self.variables[0]].load(formatter='pickle')
                ds = the_cfg.get_dataset(dataset_name)
                proc = the_cfg.get_process(ds.processes.names()[0])
               
                
                if proc.is_mc and not proc.has_tag("signal"):
                    print(f'Adding mc: {ds.processes.names()[0]}')
                    if mc_hists is None: mc_hists = the_hist
                    else: mc_hists += the_hist
                if proc.is_data:
                    print(f'Adding data: {ds.processes.names()[0]}')
                    if data_hists is None: data_hists = the_hist
                    else: data_hists += the_hist
        cfg = self.config_insts[0]
        sr_name = self.categories[0]
        aux_cats = [c for c in cfg.get_category(sr_name).x.ff_regs.values() if 'dr' in c]
        cat_list = [sr_name, *aux_cats]
        loc_list = [hist.loc(c) for c in cat_list]
        
       
        mc_h = mc_hists[{'category':loc_list, 'process' : sum, 'shift' : hist.loc('nominal')}]
        data_h = data_hists[{'category':loc_list, 'process' : sum, 'shift' : hist.loc('nominal')}]
        
        def eval_formula(formula_str, popt,make_rounding=False):
                for i,p in enumerate(popt):
                    if make_rounding:
                        formula_str = formula_str.replace(f'p{i}', '{:.3e}'.format(p))
                    else:
                        formula_str = formula_str.replace(f'p{i}',str(p))
                return formula_str
        
        #Function that performs the calculation of transfer factors
        def get_ff_corr(self, h_data, h_mc, dr_num, dr_den, name='ff_hist', label='ff_hist'):
            
            def get_single_cat(self, h, reg_name): 
                cat_name = cfg.get_category(self.categories[0]).aux['ff_regs'][reg_name]
                return h[{'category': hist.loc(cat_name)}]
            data_num = get_single_cat(self, h_data, dr_num)
            data_den = get_single_cat(self, h_data, dr_den)
            mc_num = get_single_cat(self, h_mc, dr_num)
            mc_den = get_single_cat(self, h_mc, dr_den)
            print(name)
            for nj in [0,1,2]:
                for dm in [0,1,2,10,11]:
                    print(f'DM {dm} Nj {nj}')
                    print(f"data_num: {data_num[{'tau_dm_pnet': hist.loc(dm), 'n_jets': hist.loc(nj)}].values()}")
                    print(f"data_den: {data_den[{'tau_dm_pnet': hist.loc(dm), 'n_jets': hist.loc(nj)}].values()}")
                    print(f"mc_num: {mc_num[{'tau_dm_pnet': hist.loc(dm), 'n_jets': hist.loc(nj)}].values()}")
                    print(f"mc_den: {mc_den[{'tau_dm_pnet': hist.loc(dm), 'n_jets': hist.loc(nj)}].values()}")
            num = data_num.values() - mc_num.values()

            den = data_den.values() - mc_den.values()
            ff_val = np.where((num > 0) & (den > 0),
                               num / np.maximum(den, 1),
                               -1)
            def rel_err(x):
                return x.variances()/np.maximum(x.values()**2, 1)
            
            ff_err = ff_val * ((data_num.variances() + mc_num.variances())**0.5 / np.abs(num) + (data_den.variances() + mc_den.variances())**0.5 / np.abs(den))
            
            ff_err[ff_val < 0] = 1
            h = hist.Hist.new
            for (var_name, var_axis) in cfg.x.fake_factor_method.axes.items(): 
                h = eval(f'h.{var_axis.ax_str}') 
            axes = list(h.axes[1:])
            h = h.StrCategory(['nominal', 'up', 'down'], name='syst', label='Statistical uncertainty of the fake factor')
            ff_raw = h.Weight()
            ff_raw.view().value[...,0] = ff_val
            ff_raw.view().variance[...,0] = ff_err**2
            ff_raw.name = name + '_raw'
            ff_raw.label = label + '_raw'
            
            def get_fit_model(name, dm, nj):
                if dm==0 and nj==2:
                    formula_str = 'p0+p1*x'
                    def fitf(x,p0,p1): 
                        from numpy import exp
                        return eval(formula_str)
            
                    def jac(x,p):
                        from numpy import array,exp
                        ders=array([ 1.,x])
                        return ders 
                    bounds = ([-0.5, -0.1],[0.5,0.1]) 
                else:
                    formula_str = 'p0+p1*exp(-p2*x)'
                    def fitf(x,p0,p1,p2): 
                        from numpy import exp
                        return eval(formula_str)
            
                    def jac(x,p):
                        from numpy import array,exp
                        ders=array([ 1.,
                                    exp(-p[2]*x),
                                    -1.*p[1]*x*exp(-p[2]*x)])
                        return ders
                    
                    bounds = ([-1, -5, 0],[1,5,1]) 
                return formula_str, fitf, jac, bounds
                    
        
            ff_fitted = ff_raw.copy().reset()
            ff_fitted.name = name
            ff_fitted.label = label
            
            fitres = {}
            
            dm_axis = ff_raw.axes['tau_dm_pnet']
            n_jets_axis = ff_raw.axes['n_jets']
            
            for nj in n_jets_axis:
                if nj not in fitres.keys(): fitres[nj] = {}
                for dm in dm_axis:
                    if dm not in fitres[nj].keys(): fitres[nj][dm] = {}
                    h1d = ff_raw[{'tau_dm_pnet': hist.loc(dm),
                                   'n_jets': hist.loc(nj),
                                    'syst': hist.loc('nominal')}]
                    mask = h1d.values() > 0
                    x = h1d.axes[0].centers
                    if np.sum(mask) < 2:
                        y = np.zeros_like(x)
                        y_err = np.ones_like(x)
                        x_masked = x
                    else:
                        y = h1d.values()[mask]
                        y_err = (h1d.variances()[mask])**0.5
                        x_masked = x[mask]
                        
                    formula_str, fitf, jac, the_bounds = get_fit_model(name, dm, nj)
                    popt, pcov, infodict, mesg, ier = curve_fit(fitf,
                                           x_masked,
                                           y,
                                           sigma=y_err,
                                           bounds=the_bounds,
                                           absolute_sigma=True,
                                           full_output=True
                                        )
                    print(nj,dm,(infodict['fvec'])**2)
                    fitres[nj][dm]['chi2']      = sum((infodict['fvec'])**2)
                    fitres[nj][dm]['ndf']       = len(y) - len(popt)
                    fitres[nj][dm]['popt']      = popt
                    fitres[nj][dm]['pcov']      = pcov
                    fitres[nj][dm]['x_max']     = np.max(x_masked)
                   
                    fitres[nj][dm]['jac']       = jac
                    fitres[nj][dm]['name']      = name
                    fitres[nj][dm]['fitf']      = fitf
                    fitres[nj][dm]['fitf_str']  = formula_str
                    
                    ax = h1d.axes[0]
                    binning = []
                    for idx, (x_low, x_high) in enumerate(zip(ax.edges[:-1], ax.edges[1:])):
                        if idx < ax.size-1: binning.append(np.linspace(x_low, x_high, 20, endpoint=False))
                        else: binning.append(np.linspace(x_low, x_high, 21, endpoint=True))
                    x_edges = np.concatenate(binning)
                    x_centers = (x_edges[:-1]+x_edges[1:])/2.
                    y = fitf(x_centers, *popt)
                    def err(x,jac,pcov,popt):
                            from numpy import sqrt,einsum,abs
                            return sqrt(abs(einsum('i,ij,j',jac(x,popt).T,pcov,jac(x,popt)))) 
                    import functools
                    err_y = list(map(functools.partial(err, jac=jac,pcov=pcov,popt=popt), x_centers))
                    fitres[nj][dm]['x_edges'] = x_edges
                    fitres[nj][dm]['x'] = x_centers
                    fitres[nj][dm]['y'] = {}
                    for c, shift_name in enumerate(['down', 'nom', 'up']): # if down then c=-1, if up c=+1, nominal => c=0
                        
                        fitres[nj][dm]['y'][shift_name] = np.clip(y + np.multiply(err_y,(c-1)), 0.,1.)
            return ff_raw, fitres
        
        wj_raw, wj_fitres = get_ff_corr(self,
                              data_h,
                              mc_h,
                              dr_num = 'dr_num_wj',
                              dr_den = 'dr_den_wj',
                              name='ff_wjets',
                              label='Fake factor non-QCD')
        
        qcd_raw, qcd_fitres = get_ff_corr(self,
                              data_h,
                              mc_h,
                              dr_num = 'dr_num_qcd',
                              dr_den = 'dr_den_qcd',
                              name='ff_qcd',
                              label='Fake factor QCD')
        
        
        corr_list = []            
        for fitres_per_proc in [wj_fitres, qcd_fitres]:
            binned_by_nj = []
            for nj, fitres_per_nj in fitres_per_proc.items():
                binned_by_dm = []
                for dm, fitres_per_dm in fitres_per_nj.items():
                    binned_by_syst = []
                    for the_syst, fitres in fitres_per_dm['y'].items():
                        binned_by_syst.append(cs.CategoryItem(
                            key=the_syst,
                            value=cs.Binning(
                                nodetype="binning",
                                input="tau_pt",
                                edges=fitres_per_dm['x_edges'].tolist(),
                                content=fitres.tolist(),
                                flow='clamp')
                            ))
                    binned_by_dm.append(cs.CategoryItem(
                        key=dm,
                        value=cs.Category(
                            nodetype="category",
                            input="syst",
                            content=binned_by_syst,
                            )))
                    
                binned_by_nj.append(cs.CategoryItem(
                        key=nj,
                        value=cs.Category(
                            nodetype="category",
                            input="tau_dm_pnet",
                            content=binned_by_dm,
                            )))
            corr_list.append(cs.Correction(
                name=fitres_per_proc[0][0]['name'],
                description=f"Transfer factor weights for background estimation created for mutau channel. {fitres_per_proc[0][0]['name'].split('_')[1]}",
                version=2,
                inputs=[
                    cs.Variable(name="tau_pt", type="real",description="pt of tau"),
                    cs.Variable(name="tau_dm_pnet", type="int", description="PNet decay mode of tau"),
                    cs.Variable(name="n_jets", type="int", description="Number of jets with pt > 20 GeV and eta < 4.7"),
                    cs.Variable(name="syst", type="string", description="Systematic variation 'nom', 'up', 'down'"),
                ],
                output=cs.Variable(name="weight", type="real", description="Multiplicative event weight"),
                data=cs.Category(
                    nodetype="category",
                    input="n_jets",
                    content=binned_by_nj,
                )
            ))
        cset = cs.CorrectionSet(
        schema_version=2,
        description="Fake factors",
        corrections=corr_list
        )
        self.output()['ff_json'].dump(cset.model_dump_json(exclude_unset=True), formatter="json")
        
        chi2_string = 'type nj dm chi2 ndf,'
        for fitres_per_proc in [wj_fitres, qcd_fitres]:
            for dm, fitres_per_dm in fitres_per_proc.items():
                for nj, fitres in fitres_per_dm.items():
                    chi2_string += ' '.join((fitres['name'],
                                             str(nj),
                                             str(dm),
                                             str(fitres['chi2']),
                                             str(fitres['ndf'])))
                    chi2_string += ','
        self.output()['fitres'].dump(chi2_string, formatter="json")
        
        #Plot fake factors:
        for h_name in ['wj', 'qcd']:
            h_raw       = eval(f'{h_name}_raw')
            fitres_dict = eval(f'{h_name}_fitres')
            dm_axis     = h_raw.axes['tau_dm_pnet']
            nj_axis     = h_raw.axes['n_jets']
            pt_axis     = h_raw.axes['tau_pt']
            for nj in nj_axis:
                print(f"Plotting 2d map for n jets = {nj}")
                fig, ax = plt.subplots(figsize=(12, 8))
                
                single2d_h = h_raw[{'n_jets': hist.loc(nj),
                       'syst': hist.loc('nominal')}]
                pcm = ax.pcolormesh(*np.meshgrid(*single2d_h.axes.edges), single2d_h.view().value.T, cmap="viridis", vmin=0, vmax=0.1)
                ax.set_yticks(dm_axis.centers, labels=list(map(dm_axis.bin, range(dm_axis.size))))
                plt.colorbar(pcm, ax=ax)
                plt.xlabel(single2d_h.axes.label[0])
                plt.ylabel(single2d_h.axes.label[1])
                plt.title(single2d_h.label)

                self.output()['plots']['_'.join((h_name,'nominal',f'n_jets_{str(nj)}'))].dump(fig, formatter="mpl")
                for dm in dm_axis:
                    print(f"Plotting 1d plot for n jets = {nj}, dm = {dm}")
                    h1d = h_raw[{'tau_dm_pnet': hist.loc(dm),
                                 'n_jets': hist.loc(nj),
                                    'syst': hist.loc('nominal')}]
                    fig, ax = plt.subplots(figsize=(8, 6))
                    mask = h1d.counts() > 0
                    if np.sum(mask) > 0: 
                        x = h1d.axes[0].centers[mask]
                        y = h1d.counts()[mask]
                        xerr = (np.diff(h1d.axes[0]).flatten()/2.)[mask],
                        yerr = np.sqrt(h1d.variances()).flatten()[mask],
                    else:
                        x = h1d.axes[0].centers
                        y = np.zeros_like(x)
                        xerr = (np.diff(h1d.axes[0]).flatten()/2.)
                        yerr = np.ones_like(y),
                   
                    ax.errorbar(x, y, xerr = xerr, yerr = yerr,
                                    label=f"PNet decay mode = {dm}",
                                    marker='o',
                                    fmt='o',
                                    line=None, color='#2478B7', capsize=4)
                    x = fitres_dict[nj][dm]['x']
                    y_fitf = fitres_dict[nj][dm]['y']['nom']
                    y_fitf_up = fitres_dict[nj][dm]['y']['up']
                    y_fitf_down = fitres_dict[nj][dm]['y']['down']
                    ax.plot(x,
                            y_fitf,
                            color='#FF867B')
                    ax.fill_between(x, y_fitf_up,  y_fitf_down, color='#83d55f', alpha=0.5)
                    ax.set_ylim(0,1.1*np.max(y_fitf_up))
                    #from IPython import embed; embed()
                    ax.set_xticks(pt_axis.edges, pt_axis.edges)
                    ax.set_ylabel('Fake Factor')
                    ax.set_xlabel('Tau pT [GeV]')
                    ax.grid()
                    ax.set_title(f'Jet Fake Factors : Tau PNet Decay Mode {dm}, Njets {nj}')
                    ax.annotate(rf"$\frac{{\chi^2}}{{ndf}} = \frac{{{np.round(fitres_dict[nj][dm]['chi2'],2)}}}{{{fitres_dict[nj][dm]['ndf']}}}$",
                                (0.8, 0.75),
                                xycoords='axes fraction',
                                fontsize=20)
                    
                    formula_str = eval_formula(fitres['fitf_str'],fitres_dict[nj][dm]['popt'], make_rounding=True)
                    
                    ax.annotate('y=' + formula_str,
                                (0.01, 0.95),
                                xycoords='axes fraction',
                                fontsize=12)
                    
                    self.output()['plots1d']['_'.join((h_name,str(dm),str(nj)))].dump(fig, formatter="mpl")
