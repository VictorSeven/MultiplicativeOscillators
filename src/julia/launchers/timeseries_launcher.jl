include("../core/amplitude-stochastic.jl")
include("../core/amplitude-additive.jl")
include("../core/reduced-kuramoto-cartesian.jl")

#Parameters
nharm = 30
t_thermal = 1000.0
tf = 10000.0
sys_size = 1000
w = 0.1
s2 = 0.1

if length(ARGS)==3
    #Get the program ID from Slurms's JobArray
    task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

    #  Get program parameters
    q = parse(Float64, ARGS[1])
    integration_type = ARGS[2]
    path = ARGS[3] #Where data will be stored
else
    #  Get program parameters
    q = parse(Float64, ARGS[1])
    integration_type = ARGS[2]
    path = ARGS[3] #Where data will be stored
    task_id = parse(Int, ARGS[4])
end

#Select the kind of simulation that we will do
timeseries = nothing
if integration_type=="amplitude"
    println("Amplitude selected")
    timeseries = AmplitudeEquations.get_timeseries
elseif integration_type=="additive" 
    println("Additive selected")
    timeseries = AdditiveNoise.get_timeseries
elseif  integration_type=="cartesian"
    timeseries = KuramotoCartesian.get_timeseries
end


#Run all the simulations.
#This program computes nq values from q0 to qf, and does its nrepetitions times. 
try 
    filename = "$(path)/series_$(integration_type)_$(task_id%3)"
    timeseries(nharm, t_thermal, tf, w, q, sys_size, s2, filename; nsample=100)
catch e
    println(e)
end



