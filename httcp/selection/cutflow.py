from typing import Optional
import law
from columnflow.util import maybe_import
from columnflow.selection import SelectionResult
import os

mpl = maybe_import("matplotlib")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")
np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

logger = law.logger.get_logger(__name__)

#############################################
### Plotting and saving efficiency tables ###
#############################################

def save_efficiency_table(table, dataset_name, title="Selection Yield"):
    file_path = f"/eos/user/m/mwitt/hambourg/cutflow/{title}_{dataset_name}_.txt"
    directory = os.path.dirname(file_path)
    os.makedirs(directory, exist_ok=True)

    with open(file_path, "w") as f:  # w = write
        f.write(table)
        print(f"Efficiencies saved to {file_path}")
        
def plot_eff(selections, headers, dataset_name, title="Selection Yield"): #self: Selector,
    """
    Plot efficiencies in CMS style.
    
    selections : list of lists,     Contains selection names, event counts, and efficiencies.
    headers : list of str,          Column headers where the first element is selection names.
    title : str,                    Title of the plot.
    """
    #dataset_name = self.dataset_inst.name
    mplhep.style.use("CMS")
    
    # Extract data
    step_labels = [row[0] for row in selections]  # Selection names
    event_counts = [row[1] for row in selections]  # Event counts
    efficiencies = [float(row[2]) for row in selections]  # Absolute efficiencies
    
    x_positions = np.arange(len(step_labels))
    
    # Create figure
    fig, ax1 = plt.subplots(figsize=(10.7, 10.7)) 
    
    # Plot event counts as stair plot
    edges = np.arange(len(event_counts) + 1) - 0.5
    ax1.stairs(event_counts, edges = edges , color="blue", label="Number of Events")
    ax1.set_ylabel("Event Count")
    ax1.tick_params(axis="y")
    
    # Add a second axis for efficiency
    ax2 = ax1.twinx()
    ax2.plot(x_positions, efficiencies, marker="o", color="red", linestyle="-", label="Efficiency")
    ax2.set_ylabel("Efficiency")
    ax2.tick_params(axis="y")
    ax2.set_ylim(0, 1.1)  # Efficiency should be between 0 and 1
    
    # Labeling
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(step_labels, rotation=45, ha="right")
    ax1.set_title(title, pad=50)
    
    # CMS label
    mplhep.cms.label("Private Work", loc=0)

    # Get handles and labels from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    # Combine them into one legend
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    
    # **Add efficiency values on top of each step**
    for x, eff in zip(x_positions, efficiencies):
        ax2.text(x, eff + 0.05, f"{eff:.2%}", ha="center", va="bottom", fontsize=10, color="red")


    fig.tight_layout()
    fig.savefig(f"/eos/user/m/mwitt/hambourg/cutflow/{title}_{dataset_name}.pdf", dpi=300)
    
    return fig


#############################################
###   Cutflows for event and object sels  ###
#############################################


def get_event_level_eff(events, results, dataset_name):
    from tabulate import tabulate
    steps_ = results.steps.keys()
    indiv_selections_ = []
    comb_selections_ = []
    indiv_headers_ = ["selections", "nevents", "abs eff"]
    comb_headers_ = ["selections", "nevents", "abs eff", "rel eff"]
    init = len(events)
    comb_count_array = ak.Array(np.ones(init, dtype=bool))
    for step_ in steps_:
        count_array      = results.steps[step_]
        comb_count_array = comb_count_array & count_array
        count            = ak.sum(count_array)
        comb_count       = ak.sum(comb_count_array)
        indiv_selections_.append([step_, count, round(count/init,3)])
        comb_selections_.append([step_, comb_count, round(comb_count/init,3)])
    indiv_table_ = tabulate(indiv_selections_,
                            indiv_headers_,
                            tablefmt="pretty")
    print(f"Efficiencies of individual selections: \n{indiv_table_}")
    # save the table as a .txt
    save_efficiency_table(indiv_table_, dataset_name, "Individual Selection Efficiencies")
    # plot the table as stair plot
    plot_eff(indiv_selections_, indiv_headers_, dataset_name, title="Individual Selection Efficiencies")

    comb_selections_ = np.array(comb_selections_)
    comb_selections_counts_ = comb_selections_[:,1]
    comb_den_ = np.array([init] + comb_selections_counts_[:-1].tolist())
    rel_eff_ = np.round(np.asarray(comb_selections_counts_, float)/np.asarray(comb_den_, float), decimals=3)
    comb_selections_ = np.concatenate([comb_selections_, rel_eff_[:,None]], axis=1).tolist()
    comb_table_ = tabulate(comb_selections_,
                           comb_headers_,
                           tablefmt="pretty")
    print(f"Efficiencies of combined selections: \n{comb_table_}")
    # save the table as a .txt file
    save_efficiency_table(comb_table_, dataset_name, "Combined Selection Efficiencies")
    # after tabulate inverse the y-axis order of the counts to match the other histograms
    #comb_selections_.sort(key=lambda x: x[1], reverse=True)
    # plot the table as stair plot
    plot_eff(comb_selections_, comb_headers_, dataset_name, title="Combined Selection Efficiencies")


def get_object_eff(results, tag, dataset_name):
    print(f"Object tag: {tag}")

    if tag not in results.objects:
        print(f"No object '{tag}' in SelectionResult")
        return

    obj = results.objects[tag]
    steps = obj.steps

    rows_obj = []
    rows_evt = []

    # Normalisierung: erster Schritt
    first_key = list(steps.keys())[0]
    first_mask = steps[first_key]

    n0_obj  = ak.sum(ak.sum(first_mask, axis=1))
    n0_evt  = ak.sum(ak.any(first_mask, axis=1))

    for name, mask in steps.items():
        n_obj = ak.sum(ak.sum(mask, axis=1))
        n_evt = ak.sum(ak.any(mask, axis=1))

        rows_obj.append([name, int(n_obj), round(float(n_obj / n0_obj), 3)])
        rows_evt.append([name, int(n_evt), round(float(n_evt / n0_evt), 3)])

    # Tabellen
    table_obj = tabulate(rows_obj, ["selection", f"n_{tag}", "abs eff"], tablefmt="pretty")
    table_evt = tabulate(rows_evt, ["selection", f"n_{tag}", "abs eff"], tablefmt="pretty")

    print("Object-level:\n", table_obj)
    print("Event-level:\n", table_evt)

    save_efficiency_table(table_obj, dataset_name, f"Object-Level Efficiency {tag}")
    save_efficiency_table(table_evt, dataset_name, f"Event-Level Efficiency {tag}")

    plot_eff(rows_obj, ["selection", f"n_{tag}", "abs eff"], dataset_name, title=f"Object-Level Efficiency {tag}")
    plot_eff(rows_evt, ["selection", f"n_{tag}", "abs eff"], dataset_name, title=f"Event-Level Efficiency {tag}")


def cutflow_main(events, results, dataset_name, **kwargs):
    logger.info(f"---> Inspecting event selections <---")
    get_event_level_eff(events, results, dataset_name)
    
    logger.info(f"---> Inspecting object selections <---")
    get_object_eff(results, "Muon", dataset_name)
    get_object_eff(results, "Tau", dataset_name)
    get_object_eff(results, "Jet", dataset_name)