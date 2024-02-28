include("../core/amplitude-stochastic.jl")
include("../core/amplitude-additive.jl")

#Parameters
nharm = 30
t_thermal = 100000.0
tf = 50000.0
sys_size = 1000000
s2 = 0.1

#Get the program ID from Slurms's JobArray
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

#  Get program parameters
q0 = parse(Float64, ARGS[1])
qf = parse(Float64, ARGS[2])
nq = parse(Int64,   ARGS[3])
nsims_4_diagram = parse(Int64,ARGS[4]) #Total number of launched simulations
nrepetitions = parse(Int64, ARGS[5]) 
integration_type = ARGS[6]
path = ARGS[7] #Where data will be stored

#Select the kind of simulation that we will do
phase_diagram = nothing
if integration_type=="amplitude"
    println("Amplitude selected")
    phase_diagram = AmplitudeEquations.phase_diagram
elseif integration_type=="additive" 
    println("Additive selected")
    phase_diagram = AdditiveNoise.phase_diagram
end


#Get the dq and number of simulations we have to make
n_sims_per_program = nq ÷ nsims_4_diagram
dq = (qf-q0)/nq

#The program might need, e.g. 5 parts for a phase diagram.
#This is 5 programs, each one making 100 simulations (sim0_part[0-4])
#Compute the offset so the remaining programs do sim[1-...]_part[0-4]
current_part = (task_id % nsims_4_diagram)
offset_rep = nrepetitions * (task_id ÷ nsims_4_diagram)

#Get initial and end of coupling for this program based on its task_id
q0_sim = q0 + task_id*n_sims_per_program*dq
qf_sim = q0 + ((task_id+1)*n_sims_per_program-1)*dq


#Run all the simulations.
#This program computes nq values from q0 to qf, and does its nrepetitions times. 
try 
    for index=0:nrepetitions-1
        filename = "$(path)/diagram_sim$(offset_rep + index)_part$(current_part)"
        phase_diagram(nharm, t_thermal, tf, q0_sim, qf_sim, n_sims_per_program, sys_size, s2, filename; sampling=100)
    end
catch e
    println(e)
end



