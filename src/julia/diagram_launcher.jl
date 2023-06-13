using Printf

include("core/amplitude-stochastic.jl")

using .AmplitudeEquations

#Parameters
nharm = 10
t_thermal = 500.0
tf = 100000.0
sys_size = 1e6
s2 = 0.1

task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

#  Get program parameters
q0 = parse(Float64, ARGS[1])
qf = parse(Float64, ARGS[2])
nq = parse(Int64,   ARGS[3])
nsims = parse(Int64,ARGS[4]) #Total number of launched simulations
index = parse(Int64,ARGS[5]) #Allows us to differentiate different runs 

#Get the dq and number of simulations we have to make
n_sims_per_program = nq / nsims
dq = (qf-q0)/nq

#Get initial and end of coupling for this program based on its task_id
q0_sim = q0 + task_id*n_sims_per_program*dq
qf_sim = qf + ((task_id+1)*n_sims_per_program-1)*dq

#Where data will be stored
path = "../../data/gamma"

#Run all the simulations
try 
    for (i,q) in enumerate(qs)
        filename = "$(path)/diagram_sim$(index)_part_$(i+1)"
        phase_diagram(nharm, t_thermal, tf, q0_sim, qf_sim, n_sims_per_program, sys_size, s2, filename)
    end
catch e
    println(e)
end



