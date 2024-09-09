import matplotlib.pyplot as plt
import mplhep
import pickle
import numpy as np


def create_var_plot(hists, xlabel="", ylabel="Event dencity", adaptive_ratio = None):
    """
    Create a histogram of a variable using Matplotlib with CMS style and save it to a PDF file.
    """

    # Set CMS style
    plt.style.use(mplhep.style.CMS)

    # Create Matplotlib figure and axis
    fig, (ax, ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)

    # Plotting the cutflow histogram
    color = ['red','black']
    maxy = 0
    
    for i, (name, hist) in enumerate(hists.items()):
        _ , bin_edges = hist.to_numpy()
        vals = []
        for element in hist:
            vals.append(element.count)
        local_maxy = np.max(vals/np.sum(vals))
        maxy =  local_maxy if local_maxy > maxy else maxy
        x_vals = (bin_edges[:-1] + bin_edges[1:])/2
        bin_size = (bin_edges[1] - bin_edges[0])
        #ax.scatter(x_vals, vals/np.sum(vals) , color=color[i], marker='o', alpha=0.8, label=name)
        ax.errorbar(x_vals, vals/np.sum(vals), xerr=np.ones(len(vals))*bin_size/2, fmt='o', color=color[i],alpha=0.8, label=name )
        ax.set_ylim(0,1.5*maxy)
        
        if i == 0:
            nom = np.array(vals)/np.sum(vals)
            sigma_nom = np.sqrt(np.array(vals))/np.sum(vals)
        if i == 1:
            den = np.array(vals)/np.sum(vals)
            sigma_den = np.sqrt(np.array(vals))/np.sum(vals)
            y_vals = nom/den
            y_err = (sigma_nom/nom + sigma_den/den)*y_vals
            ax_ratio.errorbar(x_vals, y_vals, yerr = y_err, xerr=np.ones(len(vals)) * bin_size / 2, fmt='o',
                                color='black', alpha=0.8, label=f"{name} / Total")
            if adaptive_ratio: ax_ratio.set_ylim(0.5*np.min(y_vals), 1.2*np.max(y_vals))
            else: ax_ratio.set_ylim(0.5, 1.5)
            ax_ratio.set_xlabel(xlabel)
            ax_ratio.set_ylabel("Ratio")
            ax_ratio.grid(True)
            ax_ratio.set_xlabel(hist.axes[0].label)
   
    # Set labels and title  
    ax.set_ylabel(ylabel)
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
    }
    mplhep.cms.label(**cms_label_kwargs)
   
    plt.tight_layout()
    
     # Save to PDF if save_path is provided
    plot_name = f'./plots/tau_joint_plot_{hist.axes[0].name}.pdf'
    fig.savefig(plot_name, bbox_inches='tight')
    print(f"Plot saved to {plot_name}")


    return fig, ax, ax_ratio

def get_hist(pickle_file):
    file_obj = open(pickle_file, 'rb')
    hist = pickle.load(file_obj)
    hist = hist.profile(axis=0)
    the_hist = hist[0,0,:] # Select only data hist withoud shift
    return the_hist

 


var_list = ["muon_pt","muon_eta","muon_phi","muon_mT","mutau_mass",
            "tau_pt","tau_eta","tau_phi",
            "jets_pt","jets_eta","jet_raw_DeepJetFlavB"]

for var in var_list:
    print(f"plotting {var}...")
    path22 = f"/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run3_2022_postEE_nano_tau_v12_limited/data_mu_f/nominal/calib__example/sel__default/prod__default/dev1/hist__{var}.pickle"
    path18 = f"/afs/cern.ch/user/s/stzakhar/work/higgs_cp/data/cf_store/analysis_higgs_cp/cf.MergeHistograms/run2_UL2018_nano_tau_v10_limited/data_ul2018_a_single_mu/nominal/calib__example/sel__default/prod__default/dev1/hist__{var}.pickle"
    
    hist18 = get_hist(path18)
    hist22 = get_hist(path22)

    create_var_plot(hists={"2022_postEE": hist22,
                           "UL2018": hist18,},)