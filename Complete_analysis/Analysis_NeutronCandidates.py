print("Iniciando el script de analisis candidatos neutrones...")
print("Importando librerias necesarias...")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import functions_spills

from collections import defaultdict

# Leer el CSV
df = pd.read_csv('/scratch/cgarcia_2002/Complete_analysis/Neutron_candidates/Second_test/neutron_candidates.csv')

# Asegurar tipos consistentes (por si event_number o start_time eran strings o enteros)
df['event_number'] = df['event_number'].astype(int)
df['start_time'] = df['start_time'].astype(float)

# Reconstruir el diccionario anidado
neutron_dict = defaultdict(lambda: defaultdict(list))

for _, row in df.iterrows():
    event_number = row['event_number']
    start_time = row['start_time']
    neutron_time = row['neutron_time']
    neutron_dict[event_number][start_time].append(neutron_time)

# Para neutron_dict
neutron_dict = {k: dict(v) for k, v in neutron_dict.items()}

df_sig = pd.read_csv('/scratch/cgarcia_2002/Complete_analysis/Neutron_candidates/Second_test/neutron_candidates_sig.csv')

df_sig['event_number'] = df_sig['event_number'].astype(int)
df_sig['start_time'] = df_sig['start_time'].astype(float)

neutron_dict_sig = defaultdict(lambda: defaultdict(list))

for _, row in df_sig.iterrows():
    event_number = row['event_number']
    start_time = row['start_time']
    neutron_time = row['neutron_time']
    neutron_dict_sig[event_number][start_time].append(neutron_time)

# Para neutron_dict_sig
neutron_dict_sig = {k: dict(v) for k, v in neutron_dict_sig.items()}

