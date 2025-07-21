print("Iniciando el script de analisis candidatos neutrones...")
print("Importando librerias necesarias...")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import functions_analysis
from scipy.optimize import curve_fit
from collections import defaultdict
import matplotlib.ticker as ticker
import glob
import os
import pickle
from tqdm import tqdm
import argparse

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['figure.figsize'] = [10, 8]
rcParams['font.size'] = 22

# Crear el parser
parser = argparse.ArgumentParser(description="Window size")

# Agregar argumento opcional con valor por defecto
parser.add_argument(
    "--window_size",
    type = int, 
    help="Window size to analyse RMS distribution",
    default=100
)

# Parsear argumentos
args = parser.parse_args()
window = args.window_size
print("Analizando window size", window)

with open('Filtered_data/datos_filtrados.pkl', 'rb') as f:
    valores_read, indices_read = pickle.load(f)

with open('Filtered_data/datos_filtrados_sig.pkl', 'rb') as f:
    valores_read_sig, indices_read_sig = pickle.load(f)

times_branch_filtered = functions_analysis.a_lista_de_arrays(valores_read, indices_read)
times_branch_filtered_sig = functions_analysis.a_lista_de_arrays(valores_read_sig, indices_read_sig)

#times_branch_filtered = times_branch_filtered[0:5000]
#times_branch_filtered_sig = times_branch_filtered_sig[0:5000]

print("Datos filtrados descargados")
N_events = len(times_branch_filtered)
N_events_sig = len(times_branch_filtered_sig)

print("Numero de eventos bkg", N_events)
print("Numero de eventos señal", N_events_sig)

print("Procesando tRMS para ventanas entre 50-300 hits")

t_RMS_list = []
for event in tqdm(np.arange(0, N_events, 1), desc="Procesando fondo"):
    threshold_times_list, _ = functions_analysis.count_nHits(times_branch_filtered[event], window, 50, 300)

    for t_in in threshold_times_list:
        t_RMS_list.append(functions_analysis.time_RMS_fun_time(times_branch_filtered[event], t_in, window))

t_RMS_list_sig = []
for event in tqdm(np.arange(0, N_events_sig, 1), desc="Procesando señal"):
    threshold_times_list, _ = functions_analysis.count_nHits(times_branch_filtered_sig[event], window, 50, 300)

    for t_in in threshold_times_list:
        t_RMS_list_sig.append(functions_analysis.time_RMS_fun_time(times_branch_filtered_sig[event], t_in, window))


"""t_RMS_list = np.concatenate([
    functions_spills.time_RMS_fun(times, window) 
    for times in tqdm(times_branch_filtered, desc="Procesando RMS time(bkg)")
])

t_RMS_list_sig = np.concatenate([
    functions_spills.time_RMS_fun(times, window) 
    for times in tqdm(times_branch_filtered_sig, desc="Procesando RMS time(sig)")
])"""

#np.savetxt(f'csv_saveData/tRMS/tRMS_{window}nsWdw_sup1.csv', t_RMS_list[~np.isnan(t_RMS_list)], delimiter=',', fmt='%d')
#np.savetxt(f'csv_saveData/tRMS/tRMS_{window}nsWdw_sig_sup1.csv', t_RMS_list_sig[~np.isnan(t_RMS_list_sig)], delimiter=',', fmt='%d')


np.savetxt(f'csv_saveData/tRMS/tRMS_{window}nsWdw_50-300.csv', t_RMS_list, delimiter=',', fmt='%d')
np.savetxt(f'csv_saveData/tRMS/tRMS_{window}nsWdw_sig_50-300.csv', t_RMS_list_sig, delimiter=',', fmt='%d')
