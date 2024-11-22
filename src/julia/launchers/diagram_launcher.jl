include("../core/amplitude-stochastic.jl")
include("../core/reduced-kuramoto-cartesian.jl")
include("../core/amplitude-additive.jl")
include("../core/kuramoto.jl")

#Parameters
nharm = 30
t_thermal = 4000.0 
tf = 15000.0 
sys_size = 100000
s2 = 0.1

nsample = 100

#Get the program ID from Slurms's JobArray
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

#  Get program parameters
q0 = parse(Float64, ARGS[1])
qf = parse(Float64, ARGS[2])
nq = parse(Int64,   ARGS[3])
nsims = parse(Int64,ARGS[4]) #Total number of launched simulations
nrepetitions = parse(Int64,ARGS[5]) 
integration_type = ARGS[6]
path = ARGS[7] #Where data will be stored

#Get the dq and number of simulations we have to make
n_sims_per_program = nq ÷ nsims
dq = (qf-q0)/nq

#Get initial and end of coupling for this program based on its task_id
q0_sim = q0 + task_id*n_sims_per_program*dq
qf_sim = q0 + ((task_id+1)*n_sims_per_program-1)*dq

#Select the kind of simulation that we will do
phase_diagram = nothing
if integration_type != "kuramoto"

    if integration_type=="amplitude"
        println("Amplitude selected")
        phase_diagram = AmplitudeEquations.phase_diagram
    elseif integration_type=="cartesian" 
        println("Cartesian selected")
        phase_diagram = KuramotoCartesian.phase_diagram
    elseif integration_type=="additive" 
        println("Additive selected")
        phase_diagram = AdditiveNoise.phase_diagram
    end

    #Run all the simulations
    try 
        for index=0:(nrepetitions-1)
            local filename = "$(path)/diagram_sim$(index)_part$(task_id)"
            phase_diagram(nharm, t_thermal, tf, q0_sim, qf_sim, n_sims_per_program, sys_size, s2, filename; sampling=nsample)
        end
    catch e
        println(e)
    end

else
    println("Kuramoto model selected")
    w = 0.1
    for index=0:(nrepetitions-1)
        local filename = "$(path)/diagram_sim$(index)_part$(task_id)"
        phase_diagram = Kuramoto.phase_diagram(sys_size, w, q0_sim, qf_sim, n_sims_per_program, s2, t_thermal, tf, filename; nsample=nsample) 
    end

end



