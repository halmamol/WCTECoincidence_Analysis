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

from matplotlib import rcParams
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams['figure.figsize'] = [10, 8]
rcParams['font.size'] = 22


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

bin_hits = 100


t_RMS_TOF = np.concatenate([
    functions_spills.time_RMS_fun(times, bin_hits) 
    for times in tqdm(times_branch_filtered, desc="Procesando RMS TOF (bkg)")
])

t_RMS_sig_TOF = np.concatenate([
    functions_spills.time_RMS_fun(times, bin_hits) 
    for times in tqdm(times_branch_filtered_sig, desc="Procesando RMS TOF (sig)")
])

np.savetxt(f'csv_saveData/TOF/TOF_100nsWdw.csv', t_RMS_TOF, delimiter=',', fmt='%d')
np.savetxt(f'csv_saveData/TOF/TOF_100nsWdw_sig.csv', t_RMS_sig_TOF, delimiter=',', fmt='%d')