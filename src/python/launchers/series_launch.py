import numpy as np
from os import system

#Simulations to be made. Offset allows to create more without overwriting previous files
nsimulations = 100 
sim_offset = 0

#Parameters for the system
kd_order = 10
s2 = 0.1
params = {"N":100, "w":0.1, "s":np.sqrt(s2), "trelax":100.0, "tf":1000.0}

#Simulate series for different q values
q_list = [0.7, 1.3]

#Make sure target folder for data exist. Classify them by N
data_path = f"../../../data/series/series_{params['N']}"
system(f"mkdir -p {data_path}")

#Compile code
system(f"g++ -O3 ../../cpp/kuramoto-mf.cpp -DMODE=SERIES -DORDER={kd_order} -o bin/kuramoto_series.out")

#Perform several simulations to average over
for sim_index in range(nsimulations):
    #For each one, do all the q_intervals
    for q in q_list:
        #Get the output path to the file and the list of arguments to the program
        output_path = f"{data_path}/series_{q:.2f}_{sim_offset + sim_index}"
        seed = 456548 + 1289*sim_index*sim_index + 12345*sim_index 
        param_string = "{N} {w} {s} {q} {trelax} {tf} {path} {seed}".format(q=q, path=output_path, **params)

        #Launch in local system or PROTEUS
        system(f"slanzarv --short --nomail bin/kuramoto_series.out {param_string}")
        #system(f"./bin/kuramoto_series.out {param_string}")