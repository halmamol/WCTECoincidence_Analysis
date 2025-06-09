import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Partition = "P3"

# Paso 1: Leer el DataFrame
df = pd.read_hdf(f"dataFrames_nHits/nHits_DataFrame_{Partition}.h5", key="df")

# Paso 2: Seleccionar la fila con event_number 1105
# Como event_number es la primera columna, podemos usar .loc o .query:
n_event = 461
fila_n= df.loc[df["event_number"] == n_event]

# fila_1105 es un DataFrame con una fila, extraemos solo las columnas con los datos
# excluimos la columna 'event_number'
datos = fila_n.drop(columns=["event_number"]).values.flatten()

# Paso 3: Graficar usando times_array (los nombres de columna)
times_array = df.columns[1:].astype(int)  # las columnas excepto 'event_number'

# Calculamos el máximo de cada fila, ignorando la columna 'event_number'
df['max_value'] = df.drop(columns=['event_number']).max(axis=1)

df_promt = df[(df["max_value"] >= 300) & (df["max_value"] <= 1000)].copy()
df_promt["idxmax"] = df_promt.drop(columns=['event_number', 'max_value']).idxmax(axis=1).astype(int)

keep_mask = []
for i, row in df_promt.iterrows():
    t_max = row["idxmax"]  # tiempo del máximo
    # Generamos los tiempos vecinos: t-2, t-1, t+1, t+2
    n=50
    vecinos = np.arange(t_max-n*1500, t_max+n*1500+1, 1500)
    vecinos = vecinos[vecinos != t_max]
    
    vecinos_validos = [v for v in vecinos if v in df_promt.drop(columns=['idxmax', 'event_number', 'max_value']).columns.astype(int)]
    vecinos_validos = [str(int(v)) for v in vecinos_validos]
    valores_vecinos = row[vecinos_validos]

    if (row[vecinos_validos] > 300).any():
        keep_mask.append(False)
    else:
        keep_mask.append(True)

# Aplicamos el filtro
df_promt = df_promt[keep_mask].reset_index(drop=True)

num_filas = df_promt.shape[0]
print("Número de filas-posibles eventos - pico entre 300 y 700 nHits y sin otros picos en 4500 ns al rededor:", num_filas)

df_promt.to_hdf(f"dataFrames_nHits/Filtered_df/nHits_DataFrame_{Partition}_Prompt.h5", key="df", mode="w")