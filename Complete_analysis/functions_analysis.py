import numpy as np
import glob
import os
import uproot

def a_lista_de_arrays(plano, indices):
    return [plano[indices[i]:indices[i+1]] for i in range(len(indices)-1)]

def count_nHits(times_branch_event_arg, window, threshold_inf, threshold_sup):
    times_branch_event = np.sort(times_branch_event_arg.copy()) #just to make sure, but it is supposed to be sorted

    if times_branch_event[0]!=times_branch_event_arg[0]:
        print("hey Carla, the intial list was not sorted : WARNING")

    threshold_times = []
    nHits_Distribution = []

    i = 0
    n = len(times_branch_event)

    while i < n:
        time_hit = times_branch_event[i]
        mask = (times_branch_event >= time_hit) & (times_branch_event < time_hit + window)

        count = mask.sum()

        if count > threshold_inf and count<threshold_sup:
            threshold_times.append(time_hit)
            nHits_Distribution.append(count)

            # Jump ahead by death window
            i += 1
            while i < n and times_branch_event[i] < time_hit + window:
                i += 1
        else:
            i += 1
            
    return threshold_times, nHits_Distribution


def get_partition_and_local_event(event_number, bordes):
    for i in range(len(bordes)-1):
        if bordes[i] <= event_number < bordes[i+1]:
            return i, event_number - bordes[i]

    print(f"Warning: Event {event_number} not found in any partition")
    return -1, -1  # Si no encuentra


if __name__ == "__main__":

    num_entries_list = []
    root_dir_bkg = "/data/cgarcia_2002/WCTE/data/2385_calib_time/"
    root_files_bkg = sorted(glob.glob(os.path.join(root_dir_bkg, "*.root")))
    root_files_bkg = sorted(root_files_bkg, key=lambda file_path: int(file_path.split("P")[-1].split(".")[0]))

    for file_path in root_files_bkg:
        print(f"Procesando archivo: {file_path}")
        file = uproot.open(file_path)
        tree = file["WCTEReadoutWindows"]
        num_entries_list.append(tree.num_entries)

    np.savetxt("Filtered_data/num_entries_list_sig.csv", num_entries_list, delimiter=",", fmt="%d")


