using Printf

include("reduced-kuramoto-cartesian.jl")

const nharm = 10
const t_thermal = 500.0
const tf = 100000.0


qs = LinRange(rc-eps, rc+eps, 100) 
const sys_size = 100
const s2 = 0.1
use_detfs = true 

sim_offset = 0
path = "../../data/series_meso/series_$(sys_size)"

task_id = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

println(task_id)
nq = length(qs) 
q = qs[task_id % nq + 1]
sim_id = task_id ÷ nq

try 
    filename = "$(path)/series_$(@sprintf "%.2f" q)_$(sim_id+sim_offset)"
    get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, filename, use_detfs)
catch e
    println(e)
end



