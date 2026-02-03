import matplotlib.pyplot as plt
import mplhep
import pickle
import numpy as np
import os

def load_histogram(pickle_file):
    """
    Load a histogram from a pickle file.
    """
    with open(pickle_file, 'rb') as f:
        hist = pickle.load(f)  # Load histogram object from file

    # Extract bin edges and values
    bin_edges = hist.axes[0].edges  # Assuming 1D histogram
    bin_values = hist.view()  # This could be multi-dimensional

    # Print the shapes of the loaded data for debugging
    print(f"Loaded histogram from {pickle_file}:")
    print(f"Bin edges shape: {bin_edges.shape}, Bin values shape: {bin_values.shape}")

    # Extract counts and ensure it's a numeric array
    if bin_values.ndim > 1:
        counts = bin_values['value'].sum(axis=tuple(range(1, bin_values.ndim)))  # Accessing the 'value' field
    else:
        counts = bin_values['value']  # In case it is already a 1D array

    # Check that the counts match the number of bins
    if counts.shape[0] != (bin_edges.shape[0] - 1):
        raise ValueError(f"Counts shape {counts.shape} does not match number of bins {bin_edges.shape[0] - 1}.")

    # Return bin edges and counts
    return bin_edges, counts


def create_var_plot(bins, counts, xlabel="", ylabel="Entries", title=""):
    """
    Create a histogram of a variable using Matplotlib with CMS style and save it to a PDF file.
    """
    plt.style.use(mplhep.style.CMS)  # Set CMS style

    # Create Matplotlib figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Calculate bin centers
    bin_centers = 0.5 * (bins[:-1] + bins[1:])  # Corrected calculation for bin centers

    # Plot the histogram
    ax.step(bin_centers, counts, where='mid', label=title, color='blue', lw=2)

    # Customize plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_title(title)
    ax.legend()
    ax.grid(True)

    # Add CMS label (default is 'Preliminary')
    mplhep.cms.label("Preliminary", lumi=8, data=True, com=13.6)

    # Save plot
    plot_name = f'./plots/{title.replace(" ", "_").lower()}_plot.pdf'
    plt.tight_layout()
    fig.savefig(plot_name, bbox_inches='tight')
    print(f"Plot saved to {plot_name}")




# Path to the histogram files
path = "/eos/user/j/jmalvaso/higgs_cp_store/analysis_httcp/cf.MergeHistograms/run3_2022_preEE_limited/dy_incl/nominal/calib__main/sel__main/prod__main/weight__main/dev/"
histogram_files = ['hist__jets_pt_no_jec.pickle', 'hist__jets_pt.pickle']

# Create full paths for the histogram files
pickle_files = [os.path.join(path, file) for file in histogram_files]

# Load histograms
bins1, counts1 = load_histogram(pickle_files[0])
bins2, counts2 = load_histogram(pickle_files[1])

# Check shapes before plotting
print(f"Shape of bins1: {bins1.shape}, counts1: {counts1.shape}")
print(f"Shape of bins2: {bins2.shape}, counts2: {counts2.shape}")

# Create the comparison plots on separate canvases
create_var_plot(bins1, counts1, xlabel="Jets pT (GeV)", ylabel="Entries", title="Jets pT No JEC")
create_var_plot(bins2, counts2, xlabel="Jets pT (GeV)", ylabel="Entries", title="Jets pT")
