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
print(f"[INFO] Analysing neutron peak partition: {Partition}")

root_file_path = f"/data/cgarcia_2002/WCTE/data/WCTE_offline_R2385S0P{Partition}.root"
# Open the ROOT file and get the TTree
file = uproot.open(root_file_path)
tree = file["WCTEReadoutWindows"] 

df_filtered = pd.read_hdf(f"dataFrames_nHits/Filtered_df/nHits_DataFrame_P{Partition}_Prompt.h5", key="df")

events_prompt = df_filtered["event_number"].tolist()
times_prompt = df_filtered["idxmax"].tolist()

i=0
Delta_t = []
for event, time in zip(events_prompt, times_prompt):
    print(f"analyzing event number {event}")
    nHits, times_array = functions_nHits.fun_window(tree, 150, event, time)
  
    # Paso 1: obtener los índices con nHits entre 15 y 100
    indices = np.where((nHits > 15) & (nHits < 100))[0]

    # Paso 2: obtener los tiempos correspondientes
    times_selected = times_array[indices]

    # Paso 3: aplicar condición adicional: time > 150
    mask = times_selected > 150

    # Paso 4: filtrar tanto los índices como los tiempos
    filtered_indices = indices[mask]
    filtered_times = times_selected[mask]

    if len(filtered_indices)!=0:

        Delta_t.append(filtered_times.max())
        print("Índices > 15:", len(filtered_indices), "in time", filtered_times)
        i+=1

print(f"Of the total number of prompt signals {len(events_prompt)}, {i} have a possible neutron peak")

bin_t = 1000
DeltaT = np.arange(0, 100, 1)

histogram_DeltaT = functions_nHits.count_nHits(np.array(Delta_t), bin_t, np.zeros((len(DeltaT))))

plt.figure()
plt.plot(DeltaT, histogram_DeltaT)
plt.xlabel(r"$\Delta t$ ($\mu$s)")
plt.ylabel("Cantidad de eventos")
plt.title(fr"Histograma de $\Delta t$ para partición {Partition}")
plt.savefig(f"histograms/DeltaT_hist_P{Partition}")