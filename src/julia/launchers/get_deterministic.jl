using DelimitedFiles

include("./core/theory_formulas.jl")
include("./core/amplitude-stochastic.jl")

using .AmplitudeEquations, .TheoryFormulas

w = 0.1
s2 = 0.1
nq = 100 
q = LinRange(0, 0.2, nq)

nsims = 10
#r = integrate_tyulkina.(w, q, s2; ntrials=nsims)
#writedlm("../../data/diagrams/theoretical/tyulkina", r)

nharm = 50 
nsims = 50
#r = integrate_full.(w, q, s2; ntrials=nsims)
#writedlm("../../data/diagrams/theoretical/deterministic_full_$nharm", r)

#nharm = 10 
#nsims = 50
#tf = 1000.0
#trelax = 500.0
#phase_diagram(nharm, trelax, tf, q[begin], q[end], nq, 1e7, s2, "../../data/diagrams/theoretical/full_system")
r = zeros(100, 2)
nsims = 100
for i=0:nsims-1
    data = readdlm("../../data/gamma/sampling_fair_10harms/diagram_sim$(i)")
    global r[:,1] += data[:,2]
    global r[:,2] += data[:,3]#data[:,2].^2
end
r /= nsims 
#r[:,2] -= r[:,1].^2
writedlm("../../data/gamma/theoretical/sf_10harms", r)

function check_angles_sto(w, q, s2; dt=0.01, tf=1000.0, nharm=30)
    angles = 2π*rand(100000)
    oldz = [mean(exp.(im*angles*k)) for k=1:nharm]
    z = Vector{ComplexF64}(undef, nharm)
    for t=0.0:dt:tf 
        full_system(nharm, w, q, s2, oldz, z, dt)
        oldz, z = z, oldz 
    end
    return angle.(oldz) 
end