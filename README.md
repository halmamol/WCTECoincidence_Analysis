Important directory Complete_Analysis. The rest were initial tests where i was only analysing one partition

1) Filtering_spills.py: save a pkl with the data filtered as datos_filtrados.pkl inside the directory Filtered_data

2) NeutronDetection.py: search for promtp + neutron candidates and save its times in the directory csv_saveData/Neutron_candidates

3) Analisis_NeutronCandidates.py: to plot the distribution of hits inside the prompt window

4) jupyter_nb/Analysis_NeutronCandidates_nb.ipynb: calculate deltaT and make the fit in the notebook reading the neutron candidates csv files

5) Analysis_PromptWindow.py: to analyse nHits distribution of windows. Create csv in csv_saveData/nHitsDistribution
    51) jupyter_nb/Analysis_PromptWindow_nb.ipynb: plots csv nHits dsitribution and ratio SIG/BKG

6) tRMS_calculation.py: calculates tRMS for windows with conditions on nHits (to calculate tRMS prompt) and saves it in a csv file in csv_saveData/tRMS
    61) jupyter_nb/tRMS_plots_nb.ipynb: plots tRMS distribution using csv files

*functions_spills.py: should include all the functions for filtering, in reality its all the functions I created in a first period of time
*functions_analysis.py: should include all the functions realted to the search of prompt-neutron candidate. In reality its all the functions I created in a second period of time

enjoy my code. sry for all the mess. XOXO 