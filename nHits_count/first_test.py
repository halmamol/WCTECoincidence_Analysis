import uproot
import numpy as np
import matplotlib.pyplot as plt

root_file_path = "/data/cgarcia_2002/WCTE/data/WCTE_offline_R2385S0P0.root"

# Open the ROOT file and get the TTree
file = uproot.open(root_file_path)
tree = file["WCTEReadoutWindows"]  

times_branch = tree["hit_pmt_times"].array(library="np")

for event in times_branch:
    div = event//1500
    res = event%1500
    for i, res_value in enumerate(res):
        if res_value != 0:
            div[i] += 1

nHits = np.zeros(179, dtype=int)

for n in div:
    nHits[int(n)] +=1

times_array = np.arange(0, 267400, 1500, dtype=int)

plt.figure()
plt.plot(times_array, nHits)
plt.show()