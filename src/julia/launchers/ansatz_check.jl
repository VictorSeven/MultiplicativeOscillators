include("../core/kuramoto.jl")

#Parameters
t_thermal = 9000.0 
tf = 10000.0 
sys_size = 100000
q = 0.3 
s2 = 0.1

filename = ARGS[1]

#Run the program
Kuramoto.get_allkd(N, q, s2, t_thermal, tf, filename)




