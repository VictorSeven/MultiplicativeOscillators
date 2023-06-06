using DelimitedFiles

include("graphs/theory_formulas.jl")
include("reduced-kuramoto-cartesian.jl")

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
r = compute_stat_r.(nharm, 1000.0, q, 100000, s2)
writedlm("../../data/diagrams/theoretical/full_stoc_$nharm", r)