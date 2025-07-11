print("Iniciando el script de detección de neutrones...")
print("Importando librerias necesarias...")

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import functions_analysis
import glob
import os
import argparse
import pickle

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['figure.figsize'] = [10, 8]
rcParams['font.size'] = 18

num_entries_list = np.loadtxt("Filtered_data/num_entries_list_bkg.csv", delimiter=",", dtype=int)
num_entries_list_sig = np.loadtxt("Filtered_data/num_entries_list_sig.csv", delimiter=",", dtype=int)

with open('Filtered_data/datos_filtrados.pkl', 'rb') as f:
    valores_read, indices_read = pickle.load(f)

with open('Filtered_data/datos_filtrados_sig.pkl', 'rb') as f:
    valores_read_sig, indices_read_sig = pickle.load(f)

times_branch_filtered = functions_analysis.a_lista_de_arrays(valores_read, indices_read)
times_branch_filtered_sig = functions_analysis.a_lista_de_arrays(valores_read_sig, indices_read_sig)

print("Datos filtrados descargados")
N_events = len(times_branch_filtered)
N_events_sig = len(times_branch_filtered_sig)

print("Numero de eventos bkg", N_events)
print("Numero de eventos señal", N_events_sig)

event_number_branch = np.arange(0, N_events, 1)
event_number_branch_sig = np.arange(0, N_events_sig, 1)
# Prompt candidates detection ###########################################################################################################

print("Buscando candidatos prompt...")
threshold_times_prompt = functions_spills.prompt_candidates(event_number_branch, times_branch_filtered, 500, 200, 100, 300)
threshold_times_prompt_sig = functions_spills.prompt_candidates(event_number_branch_sig, times_branch_filtered_sig, 500, 200, 100, 300)
print("Candidatos prompt encontrados.")

# Neutron detection ###########################################################################################################

print("Detectando neutrones en fondo...")
neutron_dict = functions_spills.neutron_detection(event_number_branch, times_branch_filtered,  threshold_times_prompt, 150000, 100, 22, threshold_sup=30, window_prompt=500)
print("Detectando neutrones en señal...")
neutron_dict_sig = functions_spills.neutron_detection(event_number_branch_sig, times_branch_filtered_sig,  threshold_times_prompt_sig, 150000, 100, 22, threshold_sup=30, window_prompt=500)

print("Prompt candidates background", sum(len(v) for v in threshold_times_prompt.values()))
print("Neutron candidates background", sum(len(v) for v in neutron_dict.values()))

print("Prompt candidates signal", sum(len(v) for v in threshold_times_prompt_sig.values()))
print("Neutron candidates signal", sum(len(v) for v in neutron_dict_sig.values()))


# Save neutron candidates to CSV files ###########################################################################################################

print("Guardando candidatos a neutrones en archivos CSV...")

bordes = np.cumsum([0] + list(num_entries_list))
bordes_sig = np.cumsum([0] + list(num_entries_list_sig))

neutron_candidates = []
for event_number, times in neutron_dict.items():
    for start_time, neutron_times in times.items():
        for neutron_time in neutron_times:
            partition, event_number_partition = functions_analysis.get_partition_and_local_event(event_number, bordes)
            neutron_candidates.append({
                'partition': partition, 
                'event_number_partition': event_number_partition,
                'event_number_total': event_number,
                'start_time': start_time,
                'neutron_time': neutron_time      
            })

neutron_candidates_sig = []
for event_number, times in neutron_dict_sig.items():
    for start_time, neutron_times in times.items():
        for neutron_time in neutron_times:
            partition, event_number_partition = functions_analysis.get_partition_and_local_event(event_number, bordes_sig)
            neutron_candidates_sig.append({
                'partition': partition, 
                'event_number_partition': event_number_partition,
                'event_number': event_number,
                'start_time': start_time,
                'neutron_time': neutron_time
            })
            
df_neutron_candidates = pd.DataFrame(neutron_candidates)
df_neutron_candidates_sig = pd.DataFrame(neutron_candidates_sig)
df_neutron_candidates.to_csv('/scratch/cgarcia_2002/Complete_analysis/Neutron_candidates/neutron_candidates_100-300_22-30_TestPartition.csv', index=False)
df_neutron_candidates_sig.to_csv('/scratch/cgarcia_2002/Complete_analysis/Neutron_candidates/neutron_candidates_sig_100-300_22-30_TestPartition.csv', index=False)

print("Archivos CSV guardados.")
print("Ejecución finalizada con éxito.")

#Plot filtered data #################################################################################################
"""
print("Generando histogramas antes y después del filtrado...")
nDetections_event_in = []
nDetections_event_fin = []

nDetections_event_in_sig = []
nDetections_event_fin_sig = []

for x in times_branch_sorted_TOF:
    nDetections_event_in.append(len(x))

for x in times_branch_filtered:
    nDetections_event_fin.append(len(x))

for x in times_branch_sorted_TOF_sig:
    nDetections_event_in_sig.append(len(x))

for x in times_branch_filtered_sig:
    nDetections_event_fin_sig.append(len(x))

n_bins = 600

hist_in, bin_edges = np.histogram(nDetections_event_in, bins=n_bins)
hist_in_sig, _ = np.histogram(nDetections_event_in_sig, bins=bin_edges)  # usa los mismos bordes

hist_filtered, _ = np.histogram(nDetections_event_fin, bins = bin_edges)
hist_filtered_sig, _ = np.histogram(nDetections_event_fin_sig, bins = bin_edges)

plt.figure(facecolor='white')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linestyle='--', linewidth=1.5, label='Background')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histogram hits/event')
plt.tight_layout()
plt.xlim(0, 4000)
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Initial_Hist_hitsEvent.png")


plt.figure(facecolor='white')
plt.fill_between(bin_edges[:-1], hist_filtered,  hatch='//////',  step='post', color='white', edgecolor='red', alpha=0.55, linestyle= '--', label='Bkg (after filter)')
plt.fill_between(bin_edges[:-1], hist_filtered_sig *N_events / N_events_sig, hatch='\\\\\\\\', step='post', color='white', edgecolor='blue', alpha=0.55, label='Data (after filter)')
plt.step(bin_edges[:-1], hist_in, where='post', color='crimson', linestyle='--', linewidth=1.5, label='Bkg (before filter)')
plt.step(bin_edges[:-1], hist_in_sig *N_events/ N_events_sig, where='post', color='navy', linestyle='-', linewidth=1.5, label='Data (before filter)')
plt.legend()
plt.xlabel('Number of Hits')
plt.ylabel('Number of Events')
plt.title('Histograms Before and After Filtering')
plt.tight_layout()
plt.xlim(0, 4000)
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Filtering.png")

bin_window = 4000

nHits_tot = functions_spills.counting_nHits_window(event_number_branch, times_branch_filtered, bin_window)
nHits_in = functions_spills.counting_nHits_window(event_number_branch, times_branch_sorted_TOF, bin_window)

nHits_tot_sig = functions_spills.counting_nHits_window(event_number_branch_sig, times_branch_filtered_sig, bin_window)
nHits_in_sig = functions_spills.counting_nHits_window(event_number_branch_sig, times_branch_sorted_TOF_sig, bin_window)

hist_in, bin_edges = np.histogram(nHits_in, bins=100, range=(0, 500))
hist_in_sig, _ = np.histogram(nHits_in_sig, bins=bin_edges)  # usa los mismos bordes

hist_filtered, _ = np.histogram(nHits_tot, bins = bin_edges)
hist_filtered_sig, _ = np.histogram(nHits_tot_sig, bins = bin_edges)

n_windows_ev = 270000 / bin_window

fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

# Top plot: Background and Signal
axs[0].step(bin_edges[:-1], hist_in / (N_events * n_windows_ev), linewidth = 1, where='post', label='Background', color='red')
axs[0].step(bin_edges[:-1], hist_in_sig / (N_events_sig * n_windows_ev), linewidth = 1, where='post', label='Data', color='blue')
axs[0].set_ylabel("Fraction of windows")
axs[0].set_xlabel(f"hits in {bin_window} ns")
axs[0].set_title("Initial data")
axs[0].set_yscale('log')
axs[0].legend()

# Bottom plot: Signal/Background Ratio
ratio = np.divide(
    hist_in_sig / (N_events_sig * n_windows_ev),
    hist_in / (N_events * n_windows_ev),
    out=np.full_like(hist_in, 0, dtype=float),
    where=hist_in > 0)
axs[1].step(bin_edges[:-1], ratio, linewidth = 1, where='post', color='green', label='Data / Background')
axs[1].set_xlabel(f"hits in {bin_window} ns")
axs[1].set_ylabel("S/B ratio")
axs[1].legend()
plt.tight_layout()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Windows_BeforeFiltering.png")

fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

# Top plot: Background and Signal
axs[0].step(bin_edges[:-1], hist_filtered / (N_events * n_windows_ev), linewidth = 1, where='post', label='Background', color='red')
axs[0].step(bin_edges[:-1], hist_filtered_sig / (N_events_sig * n_windows_ev), linewidth = 1, where='post', label='Data', color='blue')
axs[0].set_ylabel("Fraction of windows")
axs[0].set_xlabel(f"hits in {bin_window} ns")
axs[0].set_title("After Filtering total data")
axs[0].set_yscale('log')
axs[0].legend()
# Bottom plot: Signal/Background Ratio
ratio_2 = np.divide(
    hist_filtered_sig / (N_events_sig * n_windows_ev),
    hist_filtered / (N_events * n_windows_ev),
    out=np.full_like(hist_filtered, 0, dtype=float),
    where=hist_filtered > 0)
axs[1].step(bin_edges[:-1], ratio_2, linewidth = 1, where='post', color='green', label='Data / Background')
axs[1].set_xlabel(f"hits in {bin_window} ns")
axs[1].set_ylabel("S/B ratio")
axs[1].legend()
plt.tight_layout()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/Comparision_SigBkg_Windows_AfterFiltering.png")

print("Histogramas guardados.")"""




"""window_ns = 500  # podrias hacer plots con esto pero en este script no sirve de nada

times_branch_prompt = []
for event in event_number_branch:

    if event in threshold_times_prompt.keys():
    
        all_hits = [t for ref_time in threshold_times_prompt[event]
            for t in times_branch_filtered[event]
            if ref_time <= t <= ref_time + window_ns]
    else:
        all_hits = []
    
    times_branch_prompt.append(np.array(all_hits))

times_branch_prompt_sig = []
for event in event_number_branch_sig:
    if event in  threshold_times_prompt_sig.keys():

        all_hits = [t for ref_time in threshold_times_prompt_sig[event]
            for t in times_branch_filtered_sig[event]
            if ref_time <= t <= ref_time + window_ns]
    else:
        all_hits= []
    times_branch_prompt_sig.append(np.array(all_hits))"""



# Delta T calculation ###########################################################################################################

"""print("Calculando Delta T entre prompt y neutrones...")

deltaT = []
for event_number in neutron_dict:
    for start_time in neutron_dict[event_number]:
        neutron_times = neutron_dict[event_number][start_time]
        deltaT.append(min(neutron_times) - start_time)

deltaT_sig = []
for event_number in neutron_dict_sig:
    for start_time in neutron_dict_sig[event_number]:
        neutron_times = neutron_dict_sig[event_number][start_time]
        deltaT_sig.append(min(neutron_times) - start_time)


hist, bins_edges = np.histogram(deltaT, bins=100, range=(0, 100000))
hist_sig, _ = np.histogram(deltaT_sig, bins=bins_edges)

plt.figure()

plt.step(bins_edges[:-1], hist, where='post', linewidth=1, label='bkg', color='red')
plt.step(bins_edges[:-1], hist_sig*N_events/N_events_sig, where='post', linewidth=1, label='signal', color='blue')
plt.xlabel('Delta T (ns)')
plt.ylabel('Número de signal candidates')
plt.title('Delta T between neutron and prompt candidate')
plt.legend()
plt.savefig("/scratch/cgarcia_2002/Complete_analysis/Plots/DeltaT_Neutron_Prompt_100-300_22-30_TestPartitions.png")

print("Delta T calculado y gráfico guardado.")"""

