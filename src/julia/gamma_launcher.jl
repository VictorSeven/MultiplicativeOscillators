include("core/amplitude-stochastic.jl")
using .AmplitudeEquations

#include("core/reduced-kuramoto-cartesian.jl")
#using .KuramotoCartesian

#Parameters
nharm = 50
t_thermal = 500.0
tf = 10000.0
sys_size = 100000
s2 = 0.1

task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

#  Get program parameters
q0 = parse(Float64, ARGS[1])
qf = parse(Float64, ARGS[2])
nq = parse(Int64,   ARGS[3])
nsims = parse(Int64,ARGS[4]) #Total number of launched simulations
nrepetitions = parse(Int64,ARGS[5]) 
path = ARGS[6]

#Get the dq and number of simulations we have to make
n_sims_per_program = nq ÷ nsims
dq = (qf-q0)/nq

#Get initial and end of coupling for this program based on its task_id
q0_sim = q0 + task_id*n_sims_per_program*dq
qf_sim = q0 + ((task_id+1)*n_sims_per_program-1)*dq


#Run all the simulations
try 
    for index=0:nrepetitions
        filename = "$(path)/diagram_sim$(index)_part$(task_id)"
        phase_diagram(nharm, t_thermal, tf, q0_sim, qf_sim, n_sims_per_program, sys_size, s2, filename)
    end
catch e
    println(e)
end



