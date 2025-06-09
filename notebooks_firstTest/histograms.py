import uproot
import numpy as np
import matplotlib.pyplot as plt

root_file_path = "/data/cgarcia_2002/WCTE/data/WCTE_offline_R2385S0P2.root"

# Open the ROOT file and get the TTree
file = uproot.open(root_file_path)
tree = file["WCTEReadoutWindows"]  
total_events = tree.num_entries
typenames_dict = tree.typenames()

for branch_name, type_str in typenames_dict.items():
    
    if type_str == "double" or type_str=="int32_t":
        print("Creating histogram of:", branch_name)

    
        # Load the branch into a NumPy array
        values = tree[branch_name].array(library="np")  # e.g., "eventID"

        # Plot a histogram
        plt.hist(values, bins=15, edgecolor="black")
        plt.xlabel(branch_name)
        plt.ylabel("Frequency")
        #plt.grid(True)
        plt.text(0.95, 0.95, f"Events: {total_events}", 
                horizontalalignment='right', verticalalignment='top',
                transform=plt.gca().transAxes, fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
        plt.title("Partition 1")
        #plt.tight_layout()
        plt.savefig(f"Plots_histograms/P2/{branch_name}.png")