# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: /afs/cern.ch/work/a/anigamov/DesyTau/CPinHToTauTau/httcp/histogramming/main.py
# Bytecode version: 3.9.0beta5 (3425)
# Source timestamp: 2025-07-22 15:21:18 UTC (1753197678)

"""
Example event weight producer.
"""
import law
import order as od
from columnflow.histogramming import HistProducer, hist_producer
from columnflow.histogramming.default import cf_default, create_hist_from_variables, fill_hist, translate_hist_intcat_to_strcat
from columnflow.hist_util import add_hist_axis
from columnflow.columnar_util import Route
from columnflow.util import maybe_import, pattern_matcher
from columnflow.columnar_util import EMPTY_FLOAT
from columnflow.types import Any
np = maybe_import('numpy')
ak = maybe_import('awkward')
hist = maybe_import('hist')

@cf_default.hist_producer(keep_weights=None, skip_compatibility_check=True, drop_weights={"normalization_weight_inclusive"})
def httcp_hist_producer(self: HistProducer, events: ak.Array, **kwargs) -> ak.Array:
    print('Invoking httcp hist producer')
    processes = self.dataset_inst.processes.names()
    weight = ak.Array(np.ones(len(events), dtype=np.float32))
    for column in self.weight_columns:
        if ((~self.dataset_inst.has_tag("ttbar")) & ('top_pt_weight' in column)):
            print("===")
            print(weight)
            print(Route(column).apply(events),column)
            print("Skipping top_pt_weight for:", processes)
            print(weight)
            print("===")
            continue
        else :
            print("======")
            print(weight)
            weight = weight * Route(column).apply(events)
            print(column, Route(column).apply(events),column)
            print(weight)
            print("======")
    weight_dict = {}
    for cat in self.config_inst.categories.names():
        category = self.config_inst.get_category(cat)
        category_ids = events.category_ids 
        mask = ak.any(category_ids == category.id, axis = 1)
        masked_weights = ak.mask(weight, mask)
        masked_weights = ak.fill_none(masked_weights, 0.0) 
    
        masked_events = events

        if 'apply_ff' not in category.aux.keys():
            weight_dict[cat] = masked_weights
        elif category.aux['apply_ff'] == 'wj':
            print(f'applying FF weights: ff_weight_wj_nominal, category: {category.name}')
            weight_dict[cat] = masked_weights * masked_events.ff_weight_wj_nominal
        elif category.aux['apply_ff'] == 'qcd':
            print(f'applying FF weights: ff_weight_qcd_nominal, category: {category.name}')
            weight_dict[cat] = masked_weights * masked_events.ff_weight_qcd_nominal

    return (events, weight_dict)

@httcp_hist_producer.init
def httcp_hist_init(self: HistProducer) -> None:
    self.weight_columns = set()
    do_keep = pattern_matcher(self.keep_weights) if self.keep_weights else lambda _, /: True
    do_drop = pattern_matcher(self.drop_weights) if self.drop_weights else lambda _, /: False
    all_weights = self.config_inst.x.event_weights.copy()
    all_weights.update(self.dataset_inst.x('event_weights', {}))
    for weight_name, shift_insts in all_weights.items():
        if not do_keep(weight_name) or do_drop(weight_name):
            continue
        self.weight_columns.add(weight_name)
        self.uses.add(weight_name)
        self.shifts |= {shift_inst.name for shift_inst in shift_insts}

@httcp_hist_producer.create_hist
def httcp_create_hist(self: HistProducer, variables: list[od.Variable], task: law.Task, **kwargs) -> dict:
    """
    Define the histogram structure for the default histogram producer.
    Returns a dictionary of hist.Histogram objects keyed by category.
    """
    histograms = {}
    for category in self.config_inst.get_leaf_categories():
        histograms[category] = create_hist_from_variables(*variables, categorical_axes=(('process', 'intcat'), ('shift', 'intcat')), weight=True)
    return histograms

@httcp_hist_producer.fill_hist
def httcp_fill_hist(self: HistProducer, h: dict, data: dict[str, Any], task: law.Task) -> None:
    """
    Fill the histogram with the data.
    """
    print('Invoking httcp hist fill function')
    category_ids = data['category']
    data = {k: v for k, v in data.items() if k != 'category'}
    weight_dict = data['weight']
    for cat in self.config_inst.categories.names():
        category = self.config_inst.get_category(cat)
        mask = ak.any(category_ids == category.id, axis = 1)
        variable_name = [key for key in data.keys() if key != 'weight' and key != 'category' and key!= 'process' and key != 'shift']
        if len(variable_name) != 1:     
            raise ValueError(f"Expected exactly one variable to fill histogram, found {len(variable_name)}: {variable_name}")
        variable_name = variable_name[0]

        masked_variable = ak.mask(data[variable_name],mask)
        masked_variable = ak.fill_none(masked_variable, EMPTY_FLOAT)  
    
        data["weight"] = weight_dict[cat]
        data[variable_name] = masked_variable  
        fill_hist(h[category], data, last_edge_inclusive=task.last_edge_inclusive) 
     
@httcp_hist_producer.post_process_hist
def default_post_process_hist(self: HistProducer, h_dict: dict, task: law.Task) -> dict:
    """
    Post-process the his togram, converting integer to string axis for consistent lookup across configs where ids might
    be different.
    """
    hists = {}
    out_hist = None
    cat_names = [cat.name for cat in self.config_inst.get_leaf_categories()]
    for category in self.config_inst.get_leaf_categories():
        h = h_dict[category]
        axis_names = {ax.name for ax in h.axes}
        if 'process' in axis_names:
            process_map = {proc_id: self.config_inst.get_process(proc_id).name for proc_id in h.axes['process']}
            h = translate_hist_intcat_to_strcat(h, 'process', process_map)
        if 'shift' in axis_names:
            shift_map = {task.global_shift_inst.id: task.global_shift_inst.name}
            h = translate_hist_intcat_to_strcat(h, 'shift', shift_map)
        if out_hist is None:
            updated_axes = (hist.axis.StrCategory(cat_names, growth=True, name='category'),) +  h.axes 
            out_hist = hist.Hist(*updated_axes, storage=hist.storage.Weight())
            out_hist[{'category': hist.loc(category.name)}] = h.view()
        else: 
            if not h.empty():
                out_hist[{'category': hist.loc(category.name)}] = h.view()
    from IPython import embed; embed()
    return out_hist


