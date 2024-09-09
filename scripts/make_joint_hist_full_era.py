import matplotlib.pyplot as plt
import mplhep
import pickle
import numpy as np


def create_var_plot(h_dict, xlabel="", ylabel="Event dencity", adaptive_ratio = None):
    """
    Create a histogram of a variable using Matplotlib with CMS style and save it to a PDF file.
    """

    # Set CMS style
    plt.style.use(mplhep.style.CMS)

    # Create Matplotlib figure and axis
    fig, (ax, ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    # Plotting the cutflow histogram
    color = ['red','black','green','cyan']
    maxy = 0
    x_vals = h_dict["x_vals"]
    hists =  h_dict["hists"]
    bin_size = h_dict["bin_size"]
    for i, (name, hist) in enumerate(hists.items()):
        local_maxy = np.max(hist)
        maxy =  local_maxy if local_maxy > maxy else maxy
        #ax.scatter(x_vals, hist, color=color[i], marker='o', alpha=0.8, label=f"Data {name}")
        ax.errorbar(x_vals,  hist, xerr=np.ones_like(hist)*bin_size/2, fmt='o', color=color[i],alpha=0.8, label=f"Data {name}")
        ax.set_ylim(0,1.5*maxy)
        label_str = f"Events / {round(bin_size,1)}"
        if "GeV" in h_dict["label"]:
            label_str += " GeV"
        ax.set_ylabel(label_str)
        #Plot ratio histogram
        if i == 0:
            nom = hist/np.sum(hist)
            sigma_nom = np.sqrt(hist)/np.sum(hist)
        if i == 1:
            den = hist/np.sum(hist)
            sigma_den = np.sqrt(hist)/np.sum(hist)
            y_vals = nom/den
            y_err = (sigma_nom/nom + sigma_den/den)*y_vals
            ax_ratio.errorbar(x_vals, y_vals, yerr = y_err, xerr=np.ones_like(hist) * bin_size / 2, fmt='o',
                                color='black', alpha=0.8, label=f"{name} / Total")
            if adaptive_ratio: ax_ratio.set_ylim(0.5*np.min(y_vals), 1.2*np.max(y_vals))
            else: ax_ratio.set_ylim(0.8, 1.2)
            ax_ratio.set_ylabel("Ratio")
            ax_ratio.grid(True)
            ax_ratio.set_xlabel(h_dict["label"])
   
    # Set labels and title  
    # Add legend
    ax.legend()
    ax.grid(True)
    label_options = {
    "wip": "Work in progress",
    "pre": "Preliminary",
    "pw": "Private work",
    "sim": "Simulation",
    "simwip": "Simulation work in progress",
    "simpre": "Simulation preliminary",
    "simpw": "Simulation private work",
    "od": "OpenData",
    "odwip": "OpenData work in progress",
    "odpw": "OpenData private work",
    "public": "",
    }
    cms_label_kwargs = {
        "ax": ax,
        "llabel": label_options.get("pw"),
        "fontsize": 22,
        "data": True,
        "rlabel": "Data to Data"
    }
    mplhep.cms.label(**cms_label_kwargs)
   
    plt.tight_layout()
    
     # Save to PDF if save_path is provided
    plot_name = f'./plots/tau_joint_plot_{h_dict["name"]}.png'
    fig.savefig(plot_name, bbox_inches='tight')
    print(f"Plot saved to {plot_name}")


    return fig, ax, ax_ratio

def get_hist(pickle_file):
    file_obj = open(pickle_file, 'rb')
    hist = pickle.load(file_obj)
    hist = hist.profile(axis=0)
    the_hist = hist[0,0,:] # Select only data hist withoud shift
    return the_hist

def sum_up_hists(hists, groups):
    tmp_hist =  next(iter(hists.values()))
    _ , bin_edges = tmp_hist.to_numpy()
    x_vals = (bin_edges[:-1] + bin_edges[1:])/2
    bin_size = (bin_edges[1] - bin_edges[0])
    for group in groups.keys():
        summed_hist = np.zeros(len(x_vals))
        for (name, hist) in hists.items():
            if group in name:
                _ , bin_edges = hist.to_numpy()
                vals = []
                for element in hist: vals.append(element.count)
                summed_hist += np.array(vals)
        groups[group] = summed_hist
    return {"x_vals"    : x_vals,
            "hists"     : groups,
            "bin_size"  : bin_size,
            "name"      :tmp_hist.axes[0].name,
            "label"     :tmp_hist.axes[0].label}
        
 

#path22 = f"/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run3_2022_postEE_nano_tau_v12_limited/data_mu_f/nominal/calib__example/sel__default/prod__default/dev1/hist__{var}.pickle"
var_list = ["muon_pt","muon_eta","muon_phi","muon_mT","mutau_mass",
            "tau_pt","tau_eta","tau_phi"]
            #"jets_pt","jets_eta","jet_raw_DeepJetFlavB"]
eras_UL2018 = ["a","b","c","d"]
lumi_2018 = 59.83 # Inv fb

eras_2022_postEE = ["f","g"]

groups = {
    "UL2018": [],
    "2022_postEE": []
}
hists = {}
for var in var_list:
    print(f"plotting {var}...")
    for era in eras_UL2018:
        path18 = f"/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run2_UL2018_nano_tau_v10/data_ul2018_{era}_single_mu/nominal/calib__example/sel__default/prod__default/columnflow_sumbission/hist__{var}.pickle"
        tmp_hist = get_hist(path18)
        hists[f"UL2018_{era}"] = tmp_hist
    for era in eras_2022_postEE:
        path22 = f"/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run3_2022_postEE_nano_tau_v12/data_mu_{era}/nominal/calib__example/sel__default/prod__default/condor_submisson/hist__{var}.pickle"
        tmp_hist = get_hist(path18)
        hists[f"2022_postEE_{era}"] = tmp_hist
    h_dict = sum_up_hists(hists, groups)
    create_var_plot(h_dict=h_dict)