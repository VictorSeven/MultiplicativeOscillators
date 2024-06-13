using CairoMakie
using DelimitedFiles

include("style_funcs.jl")
include("palette.jl")
include("../core/theory_formulas.jl")

using .TheoryFormulas
using .ArtsyPalettes
using .StyleFuncs

"""
This auxiliar function allows one to read the simulation's data, which
is split between different files. It checks for the adequate subfolders in
`data_path`, for the simulation done with size `n`. 
The user must know how many gammas `ngammas` there are there in total per 
simulation, in how many `nparts` are they split, and the number of
simulations `nsims` for the ensemble average. 
The function overwrites the `av_r`, `av_sus` and `gammas` vectors.
"""
function read_data!(data_path, n, nparts, ngammas, nsims, av_r, av_sus, gammas; folder=nothing, ncols=4)
    #Initialize vectors
    av_r   .= zeros(ngammas) 
    av_sus .= zeros(ngammas) 
    gammas .= Vector{Float64}(undef, ngammas)

    #Read simulation data. There are several simulations
    for index=1:nsims
        data = Matrix{Float64}(undef, 0, ncols) 
        #Data is distrubuted over several parts. Put them together
        for part=1:nparts
            if isnothing(folder) 
                nextdata = readdlm("$data_path/diagrams_$(n)/diagram_$(index-1)_$(part-1)")
            else
                nextdata = readdlm("$data_path/$folder/diagram_sim$(index-1)_part$(part-1)")
            end
            data = vcat(data, nextdata) 
        end

        #Get ensemble averages over all simulations 
        av_r .+= data[1:ngammas, 2] 
        av_sus .+= data[1:ngammas, 3] 
        gammas .= data[1:ngammas, 1]
        
    end
    #Finish means
    av_r ./= nsims
    av_sus ./= nsims
end

function read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, gammas; ncols=4)
    #Initialize vectors
    av_r   .= zeros(ngammas) 
    av_sus .= zeros(ngammas) 
    gammas .= Vector{Float64}(undef, ngammas)

    #Read simulation data. There are several simulations
    for index=0:nsims-1
        data = readdlm("$data_path/diagram_sim$(index)")

        #Get ensemble averages over all simulations 
        av_r .+= data[1:ngammas, 2] 
        av_sus .+= data[1:ngammas, 3] 
        gammas .= data[1:ngammas, 1]
        
    end
    #Finish means
    av_r ./= nsims
    av_sus ./= nsims
end

function plot_thermodynamic_limit(axis, colors)
    #To read the structure of the file
    nparts = 4
    ngammas = 40 
    nsims = 100
    n = 100000
    s2 = 0.1
    data_path = "../../../data/diagrams"

    #Initialize vectors
    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas )

    #Read simulation data. There are several simulations
    read_data_simpler!("$data_path/diagrams_$n", ngammas, nsims, av_r, av_sus, gammas)
        
    #Use the theoretical formulas with a finer grid
    teogammas = LinRange(0.0, 0.2, 100)
    r   = TheoryFormulas.r_oa(teogammas, s2)
    r6  = TheoryFormulas.r_6th(teogammas, s2)
    #r2c = r_2th_cumulant(teogammas, s2)
    #r2c[teogammas .< 0.1] .= 0

    #Read solution from the numerical integration
    datatyul = readdlm("$data_path/theoretical/tyulkina")
    #data30harms; 

    #Plot finally all the stuff
    lines!(axis, teogammas, r, label="Ott-Antonsen", color=colors[1])
    lines!(axis, datatyul[:,1], datatyul[:,2], label="Tyulkina et al.", color=colors[2])
    lines!(axis, teogammas, r6, label="6th harmonic", color=colors[3])


    #Plot the simulation data
    scatter!(axis, gammas, av_r, markersize=4, color=:black, label="Simulation")


    
    #=
    nparts = 100 
    ngammas = 100 
    nsims = 100
    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas)
    #read_data!(data_path, n, nparts, ngammas, nsims, av_r, av_sus, gammas; folder="30harms_amplitude", ncols=3)
    read_data_simpler!("$data_path/30harms_amplitude", ngammas, nsims, av_r, av_sus, gammas; ncols=3)
    lines!(axis, gammas, av_sus*6, label="Full", color=:black)

    nparts = 60 
    ngammas = 60 
    nsims = 100
    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas)
    #read_data!(data_path, n, nparts, ngammas, nsims, av_r, av_sus, gammas; folder="30harms_cartesian", ncols=3)
    read_data_simpler!("$data_path/30harms_cartesian", ngammas, nsims, av_r, av_sus, gammas; ncols=3)
    lines!(axis, gammas, av_sus*6, label="Full", color=:blue)
    =#
    

    #Legend and labels
    axislegend(axis, position=(0.005, 0.9))

    axis.xlabel = L"J"
    axis.ylabel = L"R"

end

function plot_finite_sizes(axis, colormap)
    #Again, read structure of the file
    nparts = 3
    ngammas = 50
    nsims = 100
    data_path = "../../../data/diagrams"
    s2 = 0.1

    #Intialize vectors
    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas)

    #We have several system sizes. Set up colors for them
    system_sizes = [100, 1000, 10000]
    cmap =  cgrad(colormap)
    colors = [cmap.colors[i] for i in [20,50,100]]

    #Read the simulation data for each size
    i=1
    for n in system_sizes
        #read_data!(data_path, n, nparts, ngammas, nsims, av_r, av_sus, gammas)
        read_data_simpler!("$data_path/diagrams_$n", ngammas, nsims, av_r, av_sus, gammas)

        #Plot the simulation data
        scatter!(axis, gammas, av_r, markersize=4, color=colors[i], label=L"N=%$n")

        #Plot the theory line
        r_th = TheoryFormulas.finite_size_r(gammas, s2, n)
        lines!(axis, gammas, r_th, color=colors[i])
        i += 1
    end

    #Axis label and legend
    axis.xlabel = L"J"
    axislegend(axis, position=(0.005, 0.9))
end

#Start the figure with two axes 
#fig = Figure(resolution=two_col_size(2*1.618), fontsize=9, figure_padding=7)
set_theme!(StyleFuncs.two_col_figure(2*1.618))
fig = Figure(figure_padding=7)
group = fig[1,1] = GridLayout()
axs = [Axis(group[1,j], xgridvisible=false, ygridvisible=false) for j=1:2]

#Gap between plots 
#colgap!(group, 10)

#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [colors[i] for i in [2,3,1]]

#Plot in each axes
plot_thermodynamic_limit(axs[1], colors)
plot_finite_sizes(axs[2], :davos)

#Fine tuning of the decorations
hideydecorations!(axs[2])
for ax in axs
    ax.xticks = [0.0, 0.1, 0.2] 
    #ylims!(ax, 0.0, 0.85)
end

#Add axes labels
StyleFuncs.label_axes(axs)

#Save the figure
save("figure1.pdf", fig, pt_per_unit = 1)