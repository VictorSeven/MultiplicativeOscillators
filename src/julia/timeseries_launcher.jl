using Printf

include("core/reduced-kuramoto-cartesian.jl")

using .KuramotoCartesian

const nharm = 10
const t_thermal = 500.0
const tf = 10000.0

qs = [0.07, 0.1, 0.12] #0.07 for small, 0.12 for big
const sys_size = 100
const s2 = 0.1
use_detfs = true 

sim_offset = 0
path = "../../data/series_meso/series_$(sys_size)"

task_id = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 

nq = length(qs) 
q = qs[task_id % nq + 1]
sim_id = task_id ÷ nq

try 
    filename = "$(path)/series_$(@sprintf "%.2f" q)_$(sim_id+sim_offset)"
    get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, filename)
catch e
    println(e)
end



