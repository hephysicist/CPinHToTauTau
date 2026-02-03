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
import warnings
np = maybe_import('numpy')
ak = maybe_import('awkward')
hist = maybe_import('hist')

@cf_default.hist_producer(keep_weights=None,
                          skip_compatibility_check=True, 
                          drop_weights={"normalization_weight_inclusive"})
def httcp_hist_producer(self: HistProducer, events: ak.Array, **kwargs) -> ak.Array:
    print('Invoking httcp hist producer')

    processes = self.dataset_inst.processes.names()
    
    weight = ak.Array(np.ones(len(events), dtype=np.float32))

    _stitch_allow = set(self.config_inst.x.stitch_samples)

    _dataset_name = getattr(self.dataset_inst, "name", "")

    for column in self.weight_columns:
        # keep your existing skip rule for top_pt_weight
        if ((self.dataset_inst.has_tag("ttbar") ^ True) & (column == 'top_pt_weight')):
            print("===")
            print(weight)
            print(Route(column).apply(events), column)
            print("Skipping top_pt_weight for:", _dataset_name)
            print(weight)
            print("===")
            continue

        # NEW: only apply stitching_weights to the selected samples
        if column == 'stitching_weight' and _dataset_name not in _stitch_allow:
            print("===")
            print(weight)
            print(Route(column).apply(events), column)
            print("Skipping stitching_weight for:", _dataset_name)
            print(weight)
            print("===")
            continue
        # --- END NEW ---
        # default: apply the weight
        print("======")
        print(weight)
        weight = weight * Route(column).apply(events)
        print(column, Route(column).apply(events), column)
        print(weight)
        print("======")
    if 'h_ggf_htt_' in self.dataset_inst.name:
        weight = ak.where(weight<10,weight,ak.mean(weight))
    weight_dict = {}
    weight_dict['nominal'] = weight
    # weight_dict['tf_wj'] = weight * events.ff_weight_wj_nominal
    # weight_dict['tf_qcd'] = weight * events.ff_weight_qcd_nominal
    return events, weight_dict

@httcp_hist_producer.init
def httcp_hist_init(self: HistProducer) -> None:
    self.weight_columns = set()
    do_keep = pattern_matcher(self.keep_weights) if self.keep_weights else lambda _, /: True
    do_drop = pattern_matcher(self.drop_weights) if self.drop_weights else lambda _, /: False
    all_weights = self.config_inst.x.event_weights.copy()
    all_weights.update(self.dataset_inst.x('event_weights', {}))
    
    # for the_weight in self.config_inst.x.fake_factor_method.columns:
    #     self.uses.add(the_weight + '*')
        
    if self.dataset_inst.is_data: pass
    else: 
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
    for category in self.config_inst.categories.names():
        histograms[category] = create_hist_from_variables(*variables, categorical_axes=(('category', 'intcat'),('process', 'intcat'), ('shift', 'intcat')), weight=True)
    return histograms

@httcp_hist_producer.fill_hist
def httcp_fill_hist(self: HistProducer, h: dict, data: dict[str, Any], task: law.Task) -> None:
    """
    Fill the histogram with the data.
    """
    for cat_name in self.config_inst.categories.names():
        cat = self.config_inst.get_category(cat_name)
        fill_data = {}
        if 'apply_ff' not in cat.aux.keys():
            fill_data['weight'] = data['weight']['nominal']
        elif cat.aux['apply_ff'] == 'wj':
            print(f'including TF weights: ff_weight_wj_nominal, category: {cat.name}')
            fill_data['weight'] = data['weight']['tf_wj']
        elif cat.aux['apply_ff'] == 'qcd':
            print(f'applying FF weights: ff_weight_qcd_nominal, category: {cat.name}')
            fill_data['weight'] = data['weight']['tf_qcd']
        mask = ak.any(data['category'] == cat.id, axis = 1)
        fill_data['weight'] = fill_data['weight'][mask]
        fill_data['category'] = ak.full_like(fill_data['weight'], cat.id, dtype=np.int32)
        fill_data['shift'] = ak.full_like(fill_data['weight'], data['shift'], dtype=np.int32)
        fill_data['process'] = data['process'][mask]
        var_name = [v for v in data.keys() if v not in ['category','process','weight', 'shift']][0]
        fill_data[var_name] = data[var_name][mask]
        fill_hist(h[cat.name], fill_data, last_edge_inclusive=task.last_edge_inclusive) 
        print(var_name)

@httcp_hist_producer.post_process_hist
def default_post_process_hist(self: HistProducer, h_dict: dict, task: law.Task) -> dict:
    """
    Post-process the histogram, converting integer to string axis for consistent lookup across configs where ids might
    be different.
    """
    h_list = list(h_dict.values())
    h = sum(h_list[1:], h_list[0].copy())
    axis_names = {ax.name for ax in h.axes}
    if 'process' in axis_names:
        process_map = {proc_id: self.config_inst.get_process(proc_id).name for proc_id in h.axes['process']}
        h = translate_hist_intcat_to_strcat(h, 'process', process_map)
    if 'shift' in axis_names:
        shift_map = {task.global_shift_inst.id: task.global_shift_inst.name}
        h = translate_hist_intcat_to_strcat(h, 'shift', shift_map)
    if 'category' in axis_names:
        cat_map = {cat_id: self.config_inst.get_category(cat_id).name for cat_id in h.axes['category']}
        h = translate_hist_intcat_to_strcat(h, 'category', cat_map)
    return h

# ### Hist producer for fake factor estimation ###

# @cf_default.hist_producer(keep_weights=None,
#                           skip_compatibility_check=True, 
#                           drop_weights={"normalization_weight_inclusive"})
# def ff_hist_producer(self: HistProducer, events: ak.Array, **kwargs) -> ak.Array:
#     print('Invoking data-driven method hist producer')
#     weight = ak.Array(np.ones(len(events), dtype=np.float32))
    
#     for column in self.weight_columns:
#             weight = weight * Route(column).apply(events)
#     weight_dict = {}
#     weight_dict['nominal'] = weight
#     return events, weight_dict

# @ff_hist_producer.init
# def ff_hist_init(self: HistProducer) -> None:
#     self.weight_columns = set()
#     do_keep = pattern_matcher(self.keep_weights) if self.keep_weights else lambda _, /: True
#     do_drop = pattern_matcher(self.drop_weights) if self.drop_weights else lambda _, /: False
#     all_weights = self.config_inst.x.event_weights.copy()
#     all_weights.update(self.dataset_inst.x('event_weights', {}))
    
#     for the_var, var_dict in self.config_inst.x.fake_factor_method.axes.items():
#         self.uses.add(Route(var_dict.var_route))
   
#     if self.dataset_inst.is_data: pass
#     else: 
#         for weight_name, shift_insts in all_weights.items():
#             if not do_keep(weight_name) or do_drop(weight_name):
#                 continue
#             self.weight_columns.add(weight_name)
#             self.uses.add(weight_name)
#             self.shifts |= {shift_inst.name for shift_inst in shift_insts}

# @ff_hist_producer.create_hist
# def ff_create_hist(self: HistProducer, variables: list[od.Variable], task: law.Task, **kwargs) -> dict:
#     """
#     Define the histogram structure for the default histogram producer.
#     Returns a dictionary of hist.Histogram objects keyed by category.
#     """
   
#     histograms = {}
#     for category in self.config_inst.categories.names():
#         h = (hist.Hist.new.IntCat([], name="process", growth=True)
#              .IntCat([], name="category", label="category", growth=True)
#              .IntCat([], name="shift", label="shift", growth=True))
#         for (var_name, var_axis) in self.config_inst.x.fake_factor_method.axes.items(): 
#             h = eval(f'h.{var_axis.ax_str}') 
#         h = h.Weight()
#         histograms[category] = h.copy()
#     return histograms

# @ff_hist_producer.fill_hist
# def ff_fill_hist(self: HistProducer, h: dict, data: dict[str, Any], task: law.Task) -> None:
#     """
#     Fill the histogram with the data.
#     """
#     for cat_name in self.config_inst.categories.names():
#         cat = self.config_inst.get_category(cat_name)
#         fill_data = {}
#         mask = ak.any(data['category'] == cat.id, axis = 1)
#         fill_data['weight']         = data['weight']['nominal'][mask]
#         fill_data['category']       = ak.full_like(fill_data['weight'], cat.id, dtype=np.int32)
#         fill_data['shift']          = ak.full_like(fill_data['weight'], data['shift'], dtype=np.int32)
#         fill_data['process']        = data['process'][mask]
#         var_names = [v for v in data.keys() if v not in ['category','process','weight', 'shift']]
#         for the_var_name in var_names:
#             if the_var_name == 'n_jets':
#                 fill_data[the_var_name] = np.clip(data[the_var_name][mask],0,2)
#             else:
#                 fill_data[the_var_name] = data[the_var_name][mask]
#         fill_hist(h[cat.name], fill_data, last_edge_inclusive=task.last_edge_inclusive) 
     
# @ff_hist_producer.post_process_hist
# def ff_post_process_hist(self: HistProducer, h_dict: dict, task: law.Task) -> dict:
#     """
#     Post-process the his togram, converting integer to string axis for consistent lookup across configs where ids might
#     be different.
#     """
#     h_list = list(h_dict.values())
#     h = sum(h_list[1:], h_list[0].copy())
#     axis_names = {ax.name for ax in h.axes}
#     if 'process' in axis_names:
#         process_map = {proc_id: self.config_inst.get_process(proc_id).name for proc_id in h.axes['process']}
#         h = translate_hist_intcat_to_strcat(h, 'process', process_map)
#     if 'shift' in axis_names:
#         shift_map = {task.global_shift_inst.id: task.global_shift_inst.name}
#         h = translate_hist_intcat_to_strcat(h, 'shift', shift_map)
#     if 'category' in axis_names:
#         cat_map = {cat_id: self.config_inst.get_category(cat_id).name for cat_id in h.axes['category']}
#         h = translate_hist_intcat_to_strcat(h, 'category', cat_map)
#     return h
