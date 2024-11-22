include("../core/kuramoto.jl")

#Parameters
nharm = 30
w = 0.1
t_thermal = 4000.0 
tf = 15000.0 
sys_size = 100000
s2 = 0.1

#Get the program ID from Slurms's JobArray
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 
nsample = 100

#Divide our coupling in intervals with diverse number of simulations
#If the list is [q0, q1, q2 ...] intervals will be [q0, q1], [q1, q2], etc
q_intervals = [0.0, 0.07, 0.1, 0.13, 0.2]
nq_list = [10, 10, 10, 10]



#Expand the intervals
q0_list = [q for q in q_intervals[begin:end-1]]
qf_list = [q for q in q_intervals[2:end]]

#Get the dq in each one (divided with a minus one given how linrange works)
dq = @. (qf_list - q0_list) / (nq_list - 1)
#So that qf_list[end] coincides with the final point we wrote in q_intervals
dq[end] = 0

#Substract so we do not have repeated simulations when stitching intervals together
@. qf_list = qf_list - dq

#Make sure target folder for data exist. Classify them by N
data_path = "../../../data/diagrams/kuramoto_julia"

#For each one of all the q_intervals
for part in eachindex(q0_list)
    q0 = q0_list[part]
    qf = qf_list[part]
    nq = nq_list[part]

    #Get the output path to the file and the list of arguments to the program
    output_path = "$data_path/diagram_sim$(task_id)_part$part"
    phase_diagram = Kuramoto.phase_diagram(sys_size, w, q0, qf, nq, s2, t_thermal, tf, output_path; nsample=nsample) 
end


