import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_nHits 
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Partition to analyse")
parser.add_argument("partition", help="Specify the partition number") 

args = parser.parse_args()
Partition = str(args.partition)
print(f"[INFO] Analysing partition: {Partition}")

root_file_path = f"/data/cgarcia_2002/WCTE/data/WCTE_offline_R2385S0P{Partition}.root"
print(f"[INFO] Opening ROOT file: {root_file_path}")

# Open the ROOT file and get the TTree
file = uproot.open(root_file_path)
tree = file["WCTEReadoutWindows"]  

# Load branches
times_branch = tree["hit_pmt_times"].array(library="np")
event_number_branch = tree["event_number"].array(library="np")

range_hits = 1500
times_array = np.arange(0, 270000+1, range_hits, dtype=int)
nHits = np.zeros((tree.num_entries, len(times_array)))
print(f"[INFO] Initialized nHits array with shape: {nHits.shape}")

print("Processing events...")
for idx, (times, event_number) in enumerate(zip(times_branch, event_number_branch)):
    nHits[event_number, :] = functions_nHits.count_nHits(times, range_hits, nHits[event_number, :])
    if idx % 1000 == 0:
        print(f"Processed {idx+1}/{len(times_branch)} events")

print("Finished processing all events.")

# Save as DataFrame
df_nHits = pd.DataFrame(nHits, columns=[str(t) for t in times_array])
df_nHits.insert(0, "event_number", event_number_branch)

output_path = f"/home/cgarcia_2002/nHits_count/dataFrames_nHits/nHits_DataFrame_P{Partition}.h5"
df_nHits.to_hdf(output_path, key="df", mode="w")
print(f"[INFO] DataFrame saved to: {output_path}")
