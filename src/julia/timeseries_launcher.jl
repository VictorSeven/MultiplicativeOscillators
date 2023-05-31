include("reduced-kuramoto-cartesian.jl")

const nharm = 10
const t_thermal = 300.0
const tf = 1000.0

const q = 0.12 #0.07 for small, 0.12 for big
const sys_size = 1000 
const s2 = 0.1
use_detfs = true 

nthr = Threads.nthreads() 
path = "../../data/timeseries/time_mesoFS_synch_$(sys_size)"

if nthr==1
    filename = "$(path)_prueba"
    get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, filename, use_detfs)
else
    for i=1:nthr
        try 
        filename = "$(path)_$(i)"
        Threads.@spawn get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, filename, use_detfs)
        catch e
            println(e)
        end
    end
end



