using DelimitedFiles

include("core/theory_formulas.jl")
include("core/amplitude-stochastic.jl")

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
#writedlm("../../data/diagrams/theoretical/full_sys_$nharm", r)


nharm = 10 
nsims = 50
tf = 1000.0
trelax = 300.0
phase_diagram(nharm, trelax, tf, q[begin], q[end], nq, 1e7, s2, "../../data/diagrams/theoretical/full_system")

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