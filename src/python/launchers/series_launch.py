import numpy as np
from os import system

#Parameters for the system
kd_order = 1
s2 = 0.1
params = {"N":1000, "w":0.1, "s":np.sqrt(s2), "trelax":1000.0, "tf":10000.0}

#Simulate series for different q values
q_list = [0.05, 0.1, 0.2]

#Make sure target folder for data exist. Classify them by N
data_path = f"../../../data/series4dists/series_{params['N']}"
system(f"mkdir -p {data_path}")

#Compile code
system(f"g++ -O3 ../../cpp/kuramoto-mf.cpp -DMODE=SERIES -DORDER={kd_order} -o bin/kuramoto_series.out")

#Perform several simulations to average over
for q in q_list:
    print(q)
    #Get the output path to the file and the list of arguments to the program
    output_path = f"{data_path}/microscopic_{q:.2f}"
    seed = 456324548 
    param_string = "{N} {w} {s} {q} {trelax} {tf} {path} {seed}".format(q=q, seed=seed, path=output_path, **params)

    #Launch in local system or PROTEUS
    #system(f"slanzarv --short --nomail bin/kuramoto_series.out {param_string}")
    system(f"./bin/kuramoto_series.out {param_string}")