print("Iniciando el script de filtrado de spills")
print("Importando librerias necesarias")

import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills
import glob
import os
import argparse

# Crear el parser
parser = argparse.ArgumentParser(description="Partition to analyse")

# Agregar argumento opcional con valor por defecto
parser.add_argument(
    "--partition",
    help="Specify the partition number (if not specified, all partitions will be used)",
    default="all"
)

# Parsear argumentos
args = parser.parse_args()
partition = args.partition

# Mostrar resultado
if partition == "all":
    print("Analysing all partitions.")

    # Background data download ###################################################################################

    print("Cargando datos de bkg...")
    root_dir_bkg = "/data/cgarcia_2002/WCTE/data/2384_calib_time/"
    root_files_bkg = sorted(glob.glob(os.path.join(root_dir_bkg, "*.root")))

    print(f"Found {len(root_files_bkg)} background ROOT files.")

    times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, event_number_branch = functions_spills.multiple_partition(root_files_bkg)

    print("Datos de background cargados.")
    N_events = max(event_number_branch) + 1
    
    #Signal data download #############################################################################################

    print("Cargando datos de signal...")
    root_dir_sig = "/data/cgarcia_2002/WCTE/data/2385_calib_time/"
    root_files_sig = sorted(glob.glob(os.path.join(root_dir_sig, "*.root")))

    print(f"Found {len(root_files_sig)} signal ROOT files.")

    times_branch_sorted_sig, times_branch_sorted_TOF_sig, charge_branch_sorted_sig, mpmt_id_branch_sorted_sig, event_number_branch_sig = functions_spills.multiple_partition(root_files_sig)

    print("Datos de signal cargados.")
    N_events_sig = max(event_number_branch_sig) + 1

else:
    print(f"Analysing partition: {partition}")

    # Background data download ###################################################################################

    root_file_path = f"/data/cgarcia_2002/WCTE/data/2384_calib_time/WCTE_offline_R2384S0P{partition}.root"
    file = uproot.open(root_file_path)
    tree = file["WCTEReadoutWindows"]  

    times_branch_sorted, times_branch_sorted_TOF, charge_branch_sorted, mpmt_id_branch_sorted, event_number_branch = functions_spills.initial_treatment(tree)
    N_events = tree.num_entries

    #Signal data download #############################################################################################

    root_file_path = f"/data/cgarcia_2002/WCTE/data/2385_calib_time/WCTE_offline_R2385S0P{partition}.root"  #signal

    # Open the ROOT file and get the TTree
    file = uproot.open(root_file_path)
    tree_sig = file["WCTEReadoutWindows"]  

    times_branch_sorted_sig, times_branch_sorted_TOF_sig, charge_branch_sorted_sig, mpmt_id_branch_sorted_sig, event_number_branch_sig = functions_spills.initial_treatment(tree_sig)
    N_events_sig = tree_sig.num_entries

print(f"Número total de eventos en background: {N_events}")
print(f"Número total de eventos en signal: {N_events_sig}")

# Filter spills using nHits threshold ###################################################################################

print("Aplicando filtro por nHits para fondo...")
times_branch_modified, threshold_times, deleted_index_dict = functions_spills.repeat_spills_nHits(event_number_branch, times_branch_sorted_TOF, 300, 5000, 6000)
print("Aplicando filtro por nHits para señal...")
times_branch_modified_sig, threshold_times_sig, deleted_index_dict_sig = functions_spills.repeat_spills_nHits(event_number_branch_sig, times_branch_sorted_TOF_sig, 300, 5000, 6000)

charge_branch_filtered = functions_spills.delete_indices_list(charge_branch_sorted, deleted_index_dict)
charge_branch_filtered_sig = functions_spills.delete_indices_list(charge_branch_sorted_sig, deleted_index_dict_sig)
print("Filtros por nHits aplicados.")

#Filter spills using charge threshold #############################################################################################

print("Aplicando filtro por carga para fondo...")
times_branch_modified_chargesTT, charge_branch_modified_chargesTT, threshold_charges, deleted_indices = functions_spills.repeat_spills_Charge(event_number_branch, times_branch_modified, charge_branch_filtered, 100, 5000, threshold = 5000)
print("Aplicando filtro por carga para señal...")
times_branch_modified_chargesTT_sig, charge_branch_modified_chargesTT_sig, threshold_charges_sig, deleted_indices_sig = functions_spills.repeat_spills_Charge(event_number_branch_sig, times_branch_modified_sig, charge_branch_filtered_sig, 100, 5000, threshold = 5000)
print("Filtros por carga aplicados.")

# Save filtered data to pickle files #############################################################################################

import pickle

# Guardar
with open('/scratch/cgarcia_2002/Complete_analysis/Filtered_data/deleted_indices_nHits_BKG.pkl', 'wb') as f:
    pickle.dump(deleted_index_dict, f)

with open('/scratch/cgarcia_2002/Complete_analysis/Filtered_data/deleted_indices_nHits_SIG.pkl', 'wb') as f:
    pickle.dump(deleted_index_dict_sig, f)

with open('/scratch/cgarcia_2002/Complete_analysis/Filtered_data/deleted_indices_Charge_BKG.pkl', 'wb') as f:
    pickle.dump(deleted_indices, f)

with open('/scratch/cgarcia_2002/Complete_analysis/Filtered_data/deleted_indices_Charge_SIG.pkl', 'wb') as f:
    pickle.dump(deleted_indices_sig, f)
