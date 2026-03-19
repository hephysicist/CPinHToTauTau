"""
Task to produce and merge histograms.
"""

from __future__ import annotations

import itertools
import luigi
import law

from columnflow.tasks.framework.base import Requirements, AnalysisTask, DatasetTask, wrapper_factory
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
        year = self.config_inst.campaign.aux['year']
        tag = self.config_inst.campaign.aux['tag']
        channel = self.config_inst.channels.get_first().name
        return {
            "ff_json": self.target('_'.join(('fake_factors', channel, str(year), tag)) + '.json'),
            "plots": {
                f"qcd_{s}_N_b_jets_{nb}": self.target(f"fake_factor_qcd_{s}_Nbjets_{nb}.png")
                for s in ['nominal', 'up', 'down']
                for nb in [0, 1, 2]
            },
            "plots1d": {
                f"qcd_{nj}_{nb}": self.target(f"fake_factor_qcd_Njets_{nj}_Nbjets_{nb}.png")
                for nj in [0, 1, 2, 3]
                for nb in [0, 1, 2]
            },
            "fitres": self.target('_'.join(('fitres', channel, str(year), tag)) + '.json'),
        }

    @law.decorator.log
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
                print(f'Adding {ds.processes.names()[0]}')
                if proc.is_mc and not proc.has_tag("signal"):
                    if mc_hists is None: mc_hists = the_hist
                    else: mc_hists += the_hist
                if proc.is_data:  
                    if data_hists is None: data_hists = the_hist
                    else: data_hists += the_hist

        mc_h = mc_hists[{'process' : sum, 'shift' : sum}]
        data_h = data_hists[{'process' : sum, 'shift' : sum}]

        cfg = self.config_insts[0]

        def eval_formula(formula_str, popt,make_rounding=False):
                for i,p in enumerate(popt):
                    if make_rounding:
                        formula_str = formula_str.replace(f'p{i}', '{:.3e}'.format(p))
                    else:
                        formula_str = formula_str.replace(f'p{i}',str(p))
                return formula_str

        def get_ff_corr(data, mc, num_reg, den_reg, name, label):
            # Extract numerator/denominator hist
            def get_cat(h, reg):
                key = self.config_inst.get_category(self.categories[0]).aux['ff_regs'][reg]
                return h[key]

            dnum = get_cat(data, num_reg)
            dden = get_cat(data, den_reg)
            mnum = get_cat(mc, num_reg)
            mden = get_cat(mc, den_reg)

            num = dnum.values() - mnum.values()
            den = dden.values() - mden.values()
            
            ff_vals = np.where((num > 0) & (den > 0), num / np.maximum(den, 1), -1)
            ff_err = ff_vals * ((np.sqrt(dnum.variances()) + np.sqrt(mnum.variances())) / np.abs(num)
                                + (np.sqrt(dden.variances()) + np.sqrt(mden.variances())) / np.abs(den))
            ff_err[ff_vals < 0] = 1

            # Build raw histogram
            hbase = hist.Hist.new
            for ax in self.config_inst.x.fake_factor_method.axes.values():
                hbase = eval(f'hbase.{ax.ax_str}')
            hbase = hbase.StrCategory(['nominal', 'up', 'down'], name='syst', label='Stat Unc')
            raw = hbase.Weight()
            raw.view().value[..., 0] = ff_vals
            raw.view().variance[..., 0] = ff_err ** 2
            raw.name = name + '_raw'
            raw.label = label

            # Prepare fitted histogram
            fit_hist = raw.copy().reset()
            fit_hist.name = name
            fit_hist.label = label

            fit_results = {}
            for nb in raw.axes['N_b_jets']:
                fit_results[nb] = {}
                for nj in raw.axes['N_jets']:
                    slice1d = raw[{ 'N_b_jets': hist.loc(nb), 'N_jets': hist.loc(nj), 'syst': hist.loc('nominal') }]
                    x = slice1d.axes[0].centers
                    y = slice1d.values()
                    yerr = np.sqrt(slice1d.variances())
                    mask = y > 0

                    # Choose fit function
                    # if nj == 0:
                    func = lambda xx, p0, p1, p2: p0 + p1 * xx + p2 * xx ** 2
                    fstr = 'p0+p1*x+p2*x*x'
                        # bounds = ([-10, -5, -1], [10, 5, 1])
                    # else:
                    #     func = lambda xx, p0, p1, p2: p0 + p1 * np.exp(-p2 * xx)
                    #     fstr = 'p0+p1*exp(-p2*x)'
                    #     # bounds = ([-0.5, -3, 0], [0.5, 3, 0.1])

                    if mask.sum() >= 3:
                        popt, pcov = curve_fit(func, x[mask], y[mask], sigma=yerr[mask], maxfev=5000, absolute_sigma=True)
                    else:
                        popt = np.zeros(3)
                        pcov = np.zeros((3, 3))

                    # Numeric fit closure
                    fit_func = lambda xx, f=func, p=popt: f(xx, *p)
                    fit_results[nb][nj] = {'popt': popt, 'pcov': pcov, 'fstr': fstr, 'func': fit_func}

                    # Fill fitted histogram
                    for i, shift in enumerate(['down', 'nominal', 'up']):
                        vals = func(x, *(popt + (i - 1) * np.sqrt(np.diag(pcov))))
                        fit_hist.view().value[:, fit_hist.axes['N_jets'].index(nj),
                                                fit_hist.axes['N_b_jets'].index(nb),
                                                fit_hist.axes['syst'].index(shift)] = vals

            return raw, fit_hist, fit_results

        # Compute QCD fake factors
        q_raw, q_fit, q_res = get_ff_corr(data_hists, mc_hists, 'dr_num_qcd', 'dr_den_qcd', 'ff_qcd', 'Fake QCD')

        # Build and dump correction set
        corr = cs.Correction(
            name='ff_qcd', description='Fake factor QCD', version=2,
            inputs=[
                cs.Variable(name='delta_r', type='real', description='Delta R'),
                cs.Variable(name='N_jets', type='int', description='Number of jets'),
                cs.Variable(name='N_b_jets', type='int', description='Number of b jets')
            ],
            output=cs.Variable(name='weight', type='real', description='Weight'),
            data=cs.Category(
                nodetype='category', input='N_b_jets', content=[
                    cs.CategoryItem(
                        key=nb,
                        value=cs.Category(
                            nodetype='category', input='N_jets', content=[
                                cs.CategoryItem(
                                    key=nj,
                                    value=cs.Formula(
                                        nodetype='formula', variables=['delta_r'], parser='TFormula',
                                        expression=eval_formula(q_res[nb][nj]['fstr'], q_res[nb][nj]['popt'])
                                    )
                                ) for nj in q_res[nb]
                            ]
                        )
                    ) for nb in q_res
                ]
            )
        )
        cset = cs.CorrectionSet(schema_version=2, description='Fake factors', corrections=[corr])
        self.output()['ff_json'].dump(cset.json(exclude_unset=True), formatter='json')
        self.output()['fitres'].dump(str(q_res), formatter='json')

        # Plotting
        for nb in q_raw.axes['N_b_jets']:
            fig, ax = plt.subplots(figsize=(12, 8))
            h2d = q_raw[{ 'N_b_jets': hist.loc(nb), 'syst': hist.loc('nominal') }]
            pcm = ax.pcolormesh(*np.meshgrid(*h2d.axes.edges), h2d.view().value.T)
            ax.set_xlabel(h2d.axes[0].label)
            ax.set_ylabel(h2d.axes[1].label)
            plt.colorbar(pcm, ax=ax)
            self.output()['plots'][f'qcd_nominal_N_b_jets_{nb}'].dump(fig, formatter='mpl')
            for nj in q_raw.axes['N_jets']:
                fig, ax = plt.subplots(figsize=(8, 6))
                h1d = q_raw[{ 'N_jets': hist.loc(nj), 'N_b_jets': hist.loc(nb), 'syst': hist.loc('nominal') }]
                x, y = h1d.axes[0].centers, h1d.counts()
                yerr = np.sqrt(h1d.variances()).flatten()
                ax.errorbar(x, y, yerr=yerr, fmt='o', capsize=4)
                func = q_res[nb][nj]['func']
                xf = np.linspace(x.min(), x.max(), 100)
                ax.plot(xf, func(xf))
                ax.set_xlabel('Delta R')
                ax.set_ylabel('Fake Factor')
                self.output()['plots1d'][f'qcd_{nj}_{nb}'].dump(fig, formatter='mpl')