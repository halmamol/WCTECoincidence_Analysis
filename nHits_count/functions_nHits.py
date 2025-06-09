import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd   

def fun_window1(tree, bin_hits, event_number, prompt_time):

    window_tot = 100000 
    times_array = np.arange(0, window_tot+1, bin_hits, dtype=int)
    nHits = np.zeros(len(times_array))

    times_branch = tree["hit_pmt_times"].array(library="np")
    hits = times_branch[event_number]

    if (prompt_time+window_tot) >= 270000:
        try:
            hits_extra = times_branch[event_number + 1] + 270000
            print("We are extending to the next event the time range")
        except:
            hits_extra = np.array([])
            print(f"[INFO] Skipping extra window for last event ({event_number}) — no next event available.")

        hits = np.concatenate((hits, hits_extra)) 


    filtered_hits = hits[(hits > prompt_time) & (hits < prompt_time + window_tot)]

    res_in = prompt_time % bin_hits 
    if res_in==0:
        intitial_position = prompt_time // bin_hits 
    else:
        intitial_position = prompt_time // bin_hits +1 
    #intitial_position -= 1

    div = filtered_hits//bin_hits - intitial_position
    res = filtered_hits%bin_hits

    for i, res_value in enumerate(res):
        if res_value != 0:
            div[i] += 1

    for n in div:
        nHits[int(n)-1] +=1

    #print(nHits)
    for i in range(11):
        if nHits[0]<100:
            nHits = np.delete(nHits, 0)
            nHits = np.append(nHits, [0])
            """WARNING: estic posant el ultims valors buits, nomes pk totes les arrays tinguin el mateix tamany, es podria arreglar
            pero considero no necesari la diferencia de 1500ns que pot haver-hi entre uns events i altres"""
        else:
            #print("time prompt: ", prompt_time + i*bin_hits)
            return nHits, times_array
    print(f"warning, no hi ha pic a les primeres bins, check {event_number}")

    return nHits, times_array       


