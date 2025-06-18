import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  

def spill_nHitsTT(times_branch_event_arg, threshold):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    if times_branch_event[0]!=times_branch_event_arg[0]:
        print("hey Carla, the intial list was not sorted : WARNING")
    threshold_times = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + 50)
        count = mask.sum()

        if count > threshold:
            threshold_times.append(time_hit)

            # Zero out the next 50ns after the hit window
            mask_2 = (times_branch_event >= time_hit) & (times_branch_event < time_hit + 100)

            times_branch_event = times_branch_event[~mask_2]
            n = len(times_branch_event)

            # Jump ahead by 100 units
            i += 1
            while i < n and times_branch_event[i] < time_hit + 100:
                i += 1
        else:
            i += 1

    return times_branch_event, threshold_times

def spill_ChargeTT(charge_branch_event_arg, times_branch_event_arg, threshold):
    
    times_branch_event = times_branch_event_arg.copy()
    charge_branch_event = charge_branch_event_arg.copy()
    threshold_times = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + 50)
        charge_sum = charge_branch_event[mask].sum()
        
        if charge_sum > threshold:
            threshold_times.append(time_hit)

            # Zero out the next 50ns after the hit window
            mask_2 = (times_branch_event >= time_hit) & (times_branch_event < time_hit + 100)
            times_branch_event = times_branch_event[~mask_2]
            charge_branch_event = charge_branch_event[~mask_2]
            n = len(times_branch_event)

            # Jump ahead by 100 units
            i += 1
            while i < n and times_branch_event[i] < time_hit + 100:
                i += 1
        else:
            i += 1

    return times_branch_event, charge_branch_event, threshold_times

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

def plot_TotalCharge_Time(time, charge, bin_time):

    sum_charges = np.zeros(int(50/bin_time))
    div = (time-min(time))//bin_time

    for i, n in enumerate(div):
        sum_charges[int(n)] += charge[i]

    return sum_charges


def repeat_spills_nHits(event_number_branch, times_branch_sorted, threshold = 5):
    threshold_times = {}
    times_branch_modified = []

    for event in event_number_branch:
        remaining_times = times_branch_sorted[event].copy()
        all_thresholds = []

        pass_counter = 0  # Count how many times the function is applied

        # Repeat spill_nHitsTT until no more triggers are found
        while True:
            remaining_times, new_triggers = spill_nHitsTT(remaining_times, threshold)
            
            if not new_triggers:
                break

            pass_counter += 1
            all_thresholds.extend(new_triggers)

        if pass_counter > 1:
            print(f"Event {event}: spill_nHitsTT applied {pass_counter} times")


        times_branch_modified.append(remaining_times)

        if all_thresholds:
            threshold_times[event] = all_thresholds

    return times_branch_modified, threshold_times

def repeat_spills_Charge(event_number_branch, times_branch_sorted, charge_branch_sorted, threshold = 5000):
    threshold_charges = {}
    times_branch_modified_chargesTT = []
    charge_branch_modified_chargesTT = []

    for event in event_number_branch:
        times_input = times_branch_sorted[event].copy()
        charges_input = charge_branch_sorted[event].copy()

        all_thresholds = []
        pass_counter = 0

        while True:
            times_input, charges_input, threshold_charge_event = spill_ChargeTT(charges_input, times_input, threshold)

            if not threshold_charge_event:
                break

            all_thresholds.extend(threshold_charge_event)
            pass_counter += 1  

        if pass_counter > 1:
            print(f"Event {event}: spill_ChargeTT applied {pass_counter} times")

        times_branch_modified_chargesTT.append(times_input)
        charge_branch_modified_chargesTT.append(charges_input)

        if all_thresholds:
            threshold_charges[event] = all_thresholds

    return times_branch_modified_chargesTT, charge_branch_modified_chargesTT, threshold_charges