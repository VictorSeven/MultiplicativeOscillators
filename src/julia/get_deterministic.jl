using DelimitedFiles

include("graphs/theory_formulas.jl")

w = 0.1
s2 = 0.1
nq = 100 
q = LinRange(0, 0.2, nq)

ntrials = 10
r = integrate_tyulkina.(w, q, s2; ntrials=ntrials)
writedlm("../../data/diagrams/theoretical/tyulkina", r)

nharm = 50 
ntrials = 50
r = integrate_full.(w, q, s2; ntrials=ntrials)
writedlm("../../data/diagrams/theoretical/full_sys_$nharm", r)
