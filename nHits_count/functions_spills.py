import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  
import json
import awkward as ak

def spill_nHitsTT(times_branch_event_arg, threshold_inf, window, death_window, charge_branch_event = [], total = False, threshold_sup = np.inf):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    if times_branch_event[0]!=times_branch_event_arg[0]:
        print("hey Carla, the intial list was not sorted : WARNING")

    threshold_times = []
    indices_to_delete = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window)

        if len(charge_branch_event)!=0:
            count = max(charge_branch_event[mask])
        else:
            count = mask.sum()


        if count > threshold_inf and count<threshold_sup:
            threshold_times.append(time_hit)

            # Zero out the next death window ns after the hit window
            mask_2 = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window + death_window)

            indices = np.where(mask_2)[0]
            indices_to_delete.extend(indices.tolist())

            # Jump ahead by death window
            i += 1
            while i < n and times_branch_event[i] < time_hit + window + death_window:
                i += 1
        else:
            i += 1

    return threshold_times, indices_to_delete


def repeat_spills_nHits(event_number_branch, times_branch_sorted, threshold, window, death_window):
    threshold_times = {}
    times_branch_modified = []
    deleted_indices_by_event = {}

    for event in event_number_branch:
        remaining_times = times_branch_sorted[event].copy()
        all_thresholds = []
        all_deleted_indices = []
        pass_counter = 0

        while True:
            new_thresholds, indices_to_delete = spill_nHitsTT(remaining_times, threshold, window, death_window)

            if not new_thresholds:
                break

            # Save deleted indices relative to the original event array
            deleted_indices = np.where(np.isin(times_branch_sorted[event], remaining_times[indices_to_delete]))[0]
            all_deleted_indices.extend(deleted_indices.tolist())

            remaining_times = np.delete(remaining_times, indices_to_delete)
            all_thresholds.extend(new_thresholds)
            pass_counter += 1

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")

        times_branch_modified.append(remaining_times)

        if all_thresholds:
            threshold_times[event] = all_thresholds

        if all_deleted_indices:
            deleted_indices_by_event[event] = all_deleted_indices

    return times_branch_modified, threshold_times, deleted_indices_by_event


def repeat_spills_Charge(event_number_branch, times_branch_sorted, charge_branch_sorted, window, death_window, threshold = 5000):
    threshold_charges = {}
    times_branch_modified_chargesTT = []
    charge_branch_modified_chargesTT = []
    deleted_indices_by_event = {}

    for event in event_number_branch:
        remaining_times = times_branch_sorted[event].copy()
        remaining_charges = charge_branch_sorted[event].copy()
        all_thresholds = []
        all_deleted_indices = []
        pass_counter = 0

        while True:
            new_thresholds, indices_to_delete = spill_nHitsTT(remaining_times, threshold, window, death_window, remaining_charges)

            if not new_thresholds:
                break

            # Save deleted indices relative to the original event array
            deleted_indices = np.where(np.isin(times_branch_sorted[event], remaining_times[indices_to_delete]))[0]
            all_deleted_indices.extend(deleted_indices.tolist())

            remaining_times = np.delete(remaining_times, indices_to_delete)
            remaining_charges = np.delete(remaining_charges, indices_to_delete)

            all_thresholds.extend(new_thresholds)
            pass_counter += 1

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")

        times_branch_modified_chargesTT.append(remaining_times)
        charge_branch_modified_chargesTT.append(remaining_charges)


        if all_thresholds:
            threshold_charges[event] = all_thresholds

        if all_deleted_indices:
            deleted_indices_by_event[event] = all_deleted_indices


    return times_branch_modified_chargesTT, charge_branch_modified_chargesTT, threshold_charges, deleted_indices_by_event

def plot_TotalCharge_mPMT(mpmt_id, charge):

    charge_per_mpmt = {}

    for m, c in zip(mpmt_id, charge):
        if m in charge_per_mpmt:
            charge_per_mpmt[m] += c
        else:
            charge_per_mpmt[m] = c

    # Step 2: Prepare data for plotting
    mpmt_ids = sorted(charge_per_mpmt.keys())
    charges = [charge_per_mpmt[m] for m in mpmt_ids]

    # Step 3: Plot
    plt.figure(figsize=(8, 4))
    plt.bar(mpmt_ids, charges, width = 1, color='blue',  align='edge', edgecolor='navy')
    plt.xlim(0, 150)
    plt.xlabel("mPMT ID")
    plt.ylabel("Total Charge [u.a]")
    plt.title("Charge per mPMT in 50 ns window")
    plt.tight_layout()
    plt.show()

def plot_TotalCharge_Time(time, charge, bin_time, total_time):

    sum_charges = np.zeros(int(total_time/bin_time))
    div = (time-min(time))//bin_time

    for i, n in enumerate(div):
        sum_charges[int(n)] += charge[i]

    return sum_charges


def time_RMS_fun(times_event, window):

    RMS_list = []
    n = len(times_event)
    i = 0
    
    while i <n:
        t_in = times_event[i]
        mask = (times_event >= t_in) & (times_event < t_in + window)
        times_bin = times_event[mask]

        i = np.searchsorted(times_event, t_in + window, side='left')

        N = len(times_bin)
        if N == 0 or N==1:  #salen t_RMS bajos pero es mentira, no tiene ningun significado la rms de un unico valor, los evito
            RMS_list.append(np.nan)  
            continue

        mean_t = times_bin.mean()
        RMS = np.sqrt(np.mean((times_bin - mean_t) ** 2))
        RMS_list.append(RMS)

    return RMS_list


def read_mpmt_offsets(path):
    with open(path) as f:
        mcc_map = json.load(f)

        d = {}
        for k,v in zip(mcc_map.keys(), mcc_map.values()):
            card, channel = [int(i) for i in str(int(k)/100).split(".")]
            d[(card, channel)] = v

    return d

def correction_TOF(mpmt_map, mpmt_slot_branch, pmt_position, max_slot = 106, max_pos = 19):
    lookup = np.zeros((max_slot, max_pos))
    for (card, chan), shift in mpmt_map.items():
        lookup[card, chan] = shift

        # Hit Times Correction
    flat_slot_ids     = ak.ravel(mpmt_slot_branch)
    flat_pos_ids      = ak.ravel(pmt_position)
    flat_corrections  = lookup[flat_slot_ids, flat_pos_ids]
    corrections       = ak.unflatten(flat_corrections, ak.num(mpmt_slot_branch))

    return corrections

def initial_treatment(tree):

    times_branch = tree["hit_pmt_calibrated_times"].array()
    charge_branch = tree["hit_pmt_charges"].array()
    mpmt_slot_branch = tree["hit_mpmt_slot_ids"].array()
    pmt_position = tree["hit_pmt_position_ids"].array()
    event_number_branch = tree["event_number"].array(library="np")

    # Crear una máscara booleana por evento donde ambos valores sean >= 0
    valid_mask = (mpmt_slot_branch >= 0) & (pmt_position >= 0)

    # Aplicar la máscara a cada rama para eliminar los hits inválidos
    times_branch_clean = times_branch[valid_mask]
    charge_branch_clean = charge_branch[valid_mask]
    mpmt_slot_branch_clean = mpmt_slot_branch[valid_mask]
    pmt_position_clean = pmt_position[valid_mask]

    mpmt_map = read_mpmt_offsets("/home/cgarcia_2002/nHits_count/mpmt_tof_pos1.json")
    corrections = correction_TOF(mpmt_map, mpmt_slot_branch_clean, pmt_position_clean)
    corrected_times = times_branch_clean - corrections

    # Ordenar los tiempos por evento y obtener los índices
    sorted_idx = ak.argsort(corrected_times, axis=1)

    # Usar los índices para reordenar todos los branches
    times_sorted_TOF = corrected_times[sorted_idx]
    times_sorted = times_branch_clean[sorted_idx]
    charges_sorted = charge_branch_clean[sorted_idx]
    mpmt_sorted = mpmt_slot_branch_clean[sorted_idx]

    # Convertir a listas de NumPy arrays
    times_sorted_np = [np.array(evt) for evt in times_sorted]
    times_sorted_TOF_np = [np.array(evt) for evt in times_sorted_TOF]
    charges_sorted_np = [np.array(evt) for evt in charges_sorted]
    mpmt_sorted_np = [np.array(evt) for evt in mpmt_sorted]

    return times_sorted_np, times_sorted_TOF_np, charges_sorted_np, mpmt_sorted_np, event_number_branch

def delete_indices_list(list_to_delete, indices):

    new_list = []

    for event, values in enumerate(list_to_delete):
        if event in indices:
            indices_to_delete = set(indices[event])  # convert to set for faster lookup

            # Keep elements whose index is NOT in indices_to_delete
            filtered_list = [t for i, t in enumerate(values) if i not in indices_to_delete]
        else:
            filtered_list = values.copy()

        new_list.append(np.array(filtered_list))

    return new_list

def counting_nHits_window(event_number_branch, times_branch, bin_window):
    nHits = []

    for event_number in event_number_branch:
        times_branch_event = times_branch[event_number]
        n = len(times_branch_event)
        i=0

        while i<n:
            
            t_in = times_branch_event[i]
            mask = (times_branch_event >= t_in) & (times_branch_event < t_in + bin_window)
            count = mask.sum()
            nHits.append(count)
            i = np.searchsorted(times_branch_event, t_in + bin_window, side='left')

    return nHits

def prompt_candidates(event_branch, times_branch_arg, window_sliding, window_clean, threshold_inf, threshold_sup):

    def prompt_candidates_event(event, times_branch_event_arg, window_sliding, window_clean, threshold_inf, threshold_sup):
        valid_thresholds= []
        threshold_list, _ = spill_nHitsTT(times_branch_event_arg, threshold_inf, window_sliding, 0, threshold_sup = threshold_sup)

        for time_prompt in threshold_list:

            mask_clean_1 = (times_branch_event_arg >= time_prompt - window_clean) & (times_branch_event_arg < time_prompt)
            mask_clean_2 = (times_branch_event_arg > time_prompt + window_sliding) & (times_branch_event_arg < time_prompt + window_sliding + window_clean)
            
            if mask_clean_1.sum() != 0 :
                clean_list, _ = spill_nHitsTT(times_branch_event_arg[mask_clean_1], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup)
            else:
                clean_list = []
            
            if mask_clean_2.sum() != 0 :
                clean_list_2, _ = spill_nHitsTT(times_branch_event_arg[mask_clean_2], threshold_inf, window_sliding, 0, threshold_sup=threshold_sup)
            else:
                clean_list_2 = []

            if len(clean_list) + len(clean_list_2) == 0:
                valid_thresholds.append(time_prompt)
            else:
                
                print(f"Not using trigger {time_prompt} because it has other signal_candidate too close in event {event}") 

        return valid_thresholds
    
    theshold_times = {}
    
    for event in event_branch:
        threshold_times_event = prompt_candidates_event(event, times_branch_arg[event], window_sliding, window_clean, threshold_inf, threshold_sup)
        
        if len(threshold_times_event)!=0:
            theshold_times[event] = threshold_times_event

    return theshold_times