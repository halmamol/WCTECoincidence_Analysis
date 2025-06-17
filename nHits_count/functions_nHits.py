import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd   

def fun_window(tree, bin_hits, event_number, prompt_time):

    window_tot = 100000 
    times_array = np.arange(0, window_tot+1, bin_hits, dtype=int)
    nHits = np.zeros(len(times_array))

    times_branch = tree["hit_pmt_calibrated_times"].array(library="np")
    hits = times_branch[event_number]

    if (prompt_time+window_tot) >= 270000:
        try:
            hits_extra = times_branch[event_number + 1] + 270000
            print("We are extending to the next event the time range")
        except:
            hits_extra = np.array([])
            print(f"[INFO] Skipping extra window for last event ({event_number}) — no next event available.")

        hits = np.concatenate((hits, hits_extra)) 

    filtered_hits = hits[(hits >= prompt_time) & (hits <= prompt_time + window_tot)]

    adjusted_counts = filtered_hits - prompt_time #cuidado prompt_time(1500) tiene que ser multiple de bins? sino calcular bin y restarle? dudas

    # Luego llamas a la función con esos counts ajustados
    nHits = count_nHits(adjusted_counts, bin_hits, nHits)
    n = 0

    for i in range(int(1500/bin_hits)):
        
        if nHits[0]<(0.2*bin_hits): #esto se tiene ajustar al bin tmb
            nHits = np.delete(nHits, 0)
            nHits = np.append(nHits, [0])
            n += 1
        else:
            if n != 0:
                print(f"[INFO] {n} bins removed from the beginning of the histogram for event {event_number}")
                new_hits = hits[(hits > prompt_time + window_tot) & (hits <= prompt_time + window_tot + n*bin_hits)]
                adjusted_new_counts = new_hits - (prompt_time + n*bin_hits)
                nHits = count_nHits(adjusted_new_counts, bin_hits, nHits)
                
            #print("time prompt: ", prompt_time + i*bin_hits)
            return nHits, times_array
    print(f"warning, no hi ha pic a les primeres bins, check {event_number}")

    return nHits, times_array       

def count_nHits(counts, bin_hist, histogram):

    div = counts//bin_hist 
  
    for n in div:
        histogram[int(n)] +=1

    return histogram.astype(int)

def filter_neighbor(df, n):
    """aqui filtramos para que no haya ningun pico cercano alto"""
    print(f"[INFO] Número de eventos antes del filtro: {len(df)}")
    
    df = df.reset_index(drop=True)
    keep_mask = np.ones(len(df), dtype=bool)  # asumimos todos True inicialmente
    cols = df.drop(columns=['idxmax', 'event_number', 'max_value']).columns
    cols_int = set(cols.astype(int))
    

    for i, row in df.iterrows():
        t_max = row["idxmax"]  # tiempo del máximo

        # Generamos los tiempos vecinos: t-2, t-1, t+1, t+2
        vecinos = np.arange(t_max-n*1500, t_max+n*1500+1, 1500)
        vecinos = vecinos[vecinos != t_max]

        vecinos_validos = [str(int(v)) for v in vecinos if v in cols_int]
        
        if (row[vecinos_validos] > 300).any():
            print(f"Evento {int(row['event_number'])} descartado por vecinos con valor > 300.")
            keep_mask[i] = False
            continue 

        if t_max+n*1500>270000 and i + 1 < len(df):
            print(f"Evento {int(row['event_number'])}, accediendo valores evento {int(row['event_number'])+1}")
            up_vecinos = np.arange(0, (t_max+n*1500 - 270000)+1, 1500)
            vecinos_validos_up = [str(int(v)) for v in up_vecinos if v in cols_int]
            if (df.loc[i + 1, vecinos_validos_up] > 300).any():
                print(f"Evento {int(row['event_number'])} descartado por vecinos superior con valor > 300.")
                keep_mask[i] = False
                continue

        elif t_max-n*1500<0 and i - 1 >= 0:
            print(f"Evento {int(row['event_number'])}, accediendo valores evento {int(row['event_number'])-1}")
            down_vecinos = np.arange(270000 + t_max-n*1500, 27000+1, 1500)
            vecinos_validos_down = [str(int(v)) for v in down_vecinos if v in cols_int]
            if (df.loc[i -1, vecinos_validos_down] > 300).any():
                print(f"Evento {int(row['event_number'])} descartado por vecinos inferior con valor > 300.")
                keep_mask[i] = False
                continue
 

    # Aplicamos el filtro
    df_filtrado = df[keep_mask].reset_index(drop=True)
    print(f"[INFO] Número de eventos después del filtro: {len(df_filtrado)}")
    return df_filtrado