import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_nHits 
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Partition to analyse")
parser.add_argument("partition", help="Specify the partition number") 

args = parser.parse_args()
Partition = str(args.partition)
print(f"[INFO] Filtering partition: {Partition}")

#Download dataframe
df = pd.read_hdf(f"dataFrames_nHits/nHits_DataFrame_P{Partition}.h5", key="df")

#Buscamos punto maximo del histograma nHits
df['max_value'] = df.drop(columns=['event_number']).max(axis=1)

total_eventos = df.shape[0]
count_0_300 = df[(df["max_value"] < 300)].shape[0]
count_300_700 = df[(df["max_value"] >= 300) & (df["max_value"] <= 700)].shape[0]
count_700_up = df[df["max_value"] > 700].shape[0]

# Mostrar resultados
print(f"Máximo entre 0–299:     {count_0_300}/{total_eventos} eventos")
print(f"Máximo entre 300–700:  {count_300_700}/{total_eventos} eventos")
print(f"Máximo > 701:          {count_700_up}/{total_eventos} eventos")

#llamamos df prompt a los eventos que tienen picos entre 300 y 700 hits
df_prompt = df[(df["max_value"] >= 300) & (df["max_value"] <= 700)].copy()

#idxmax es el tiempo en el que el pico maximo se ha dado
df_prompt["idxmax"] = df_prompt.drop(columns=['event_number', 'max_value']).idxmax(axis=1).astype(int)

#filtro para que el pico sea limpio - no hay picos en 15000ns cercanos
df_prompt = functions_nHits.filter_neighbor(df_prompt, 50)
num_filas = df_prompt.shape[0]

print("Número de posibles eventos - pico entre 300 y 700 nHits y sin otros picos en 15000 ns al rededor:", num_filas)

df_prompt.to_hdf(f"dataFrames_nHits/Filtered_df/nHits_DataFrame_P{Partition}_Prompt.h5", key="df", mode="w")