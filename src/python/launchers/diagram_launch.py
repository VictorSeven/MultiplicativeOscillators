import numpy as np
from os import system

#Simulations to be made. Offset allows to create more without overwriting previous files
nsimulations = 100
sim_offset = 0

#Parameters for the system
s2 = 0.1 #Location of critical point
kd_order = 1
params = {"N":100, "w":0.1, "s":np.sqrt(s2), "trelax":100.0, "tf":1000.0}

#Divide our coupling in intervals with diverse number of simulations
#If the list is [q0, q1, q2 ...] intervals will be [q0, q1], [q1, q2], etc
q_intervals = [0.0, 0.07, 0.13, 0.20] 
nq_list = [10, 30, 10]

#Expand the intervals
q0_list = [q for q in q_intervals[:-1]]
qf_list = [q for q in q_intervals[1:]]

#Make sure target folder for data exist. Classify them by N
data_path = f"../../../data/diagrams/diagrams_{params['N']}"
system(f"mkdir -p {data_path}")

#Compile code
system(f"g++ -O3 ../../cpp/kuramoto-mf.cpp -DMODE=DIAGRAM -DORDER={kd_order} -o bin/kuramoto_diagram.out")

#Perform several simulations to average over
for sim_index in range(nsimulations):
    #For each one, do all the q_intervals
    for part,(q0, qf, nq) in enumerate(zip(q0_list, qf_list, nq_list)):
        #Get the output path to the file and the list of arguments to the program
        output_path = f"{data_path}/diagram_{sim_offset + sim_index}_{part}"
        seed = 456548 + 1289*sim_index*sim_index + 12345*sim_index 
        param_string = "{N} {w} {s} {q0} {qf} {nq} {trelax} {tf} {path} {seed}".format(q0=q0, qf=qf, nq=nq, path=output_path, seed=seed, **params)

        #Launch in local system or PROTEUS
        system(f"slanzarv --short --nomail bin/kuramoto_diagram.out {param_string}")
        #system(f"./bin/kuramoto_diagram.out {param_string}")