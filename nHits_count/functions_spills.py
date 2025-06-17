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
            mask_2 = (times_branch_event >= time_hit + 50) & (times_branch_event < time_hit + 100)

            #if np.any(mask_2):
            #    print(f"Removing {mask_2.sum()} elements between {time_hit + 50} and {time_hit + 100} ns at ={i}")
    
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
            mask_2 = (times_branch_event >= time_hit + 50) & (times_branch_event < time_hit + 100)
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
    # Step 1: Sum charges per mPMT (manual dict logic)
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
    plt.bar(mpmt_ids, charges, width = 0.5, color='green',  align='edge', edgecolor='navy')
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

