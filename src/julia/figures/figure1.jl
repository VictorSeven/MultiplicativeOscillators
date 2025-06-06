using CairoMakie
using Colors
using DelimitedFiles

include("style_funcs.jl")
include("palette.jl")
include("../core/theory_formulas.jl")

using .TheoryFormulas
using .ArtsyPalettes
using .StyleFuncs


"""
Reads the data at data_path to obtain mean and variance stored in the files
"""
function read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, gammas)
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

function plot_thermodynamic_limit(axis, colors, s2, data_folder; ngammas=43, nsims=100, eqnumber="Eq. (9)")
    #To read the structure of the file
    data_path = "../../../data/diagrams"

    #Initialize vectors
    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas )

    #Read simulation data. There are several simulations
    #read_data_simpler!("$data_path/diagrams_$n", ngammas, nsims, av_r, av_sus, gammas)
    read_data_simpler!("$data_path/$data_folder", ngammas, nsims, av_r, av_sus, gammas)
        
    #Use the theoretical formulas with a finer grid
    teogammas = LinRange(0.0, 2*s2, 100)
    r   = TheoryFormulas.r_oa(teogammas, s2)
    r6  = TheoryFormulas.r_6th(teogammas, s2)
    rtyul = TheoryFormulas.diagram_tyulkina(teogammas, 0.1, s2)

    #Plot finally all the stuff
    lines!(axis, teogammas, r, label="Ott-Antonsen", color=colors[1])
    lines!(axis, teogammas, rtyul, label="Tyulkina et al.", color=colors[2])
    lines!(axis, teogammas, r6, label=eqnumber, color=colors[3])


    #Plot the simulation data (for paper I set markersize manually to 4)
    scatter!(axis, gammas, av_r, color=:black, label="Simulation")


    #Legend and labels
    axislegend(axis, position=(0.005, 0.95))

    axis.xlabel = L"J"
    axis.ylabel = L"R"

end

function plot_finite_sizes(axis, basecolor; use_OA=false)
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

    #Create a nice shading from the selected color's hue
    fshue = LCHuv(parse(Colorant, basecolor)).h
    fscolors = sequential_palette(fshue, 100; b=0.9, s=0.7)
    fscolors = [fscolors[i] for i in [60, 80, 100]]

    #Create a nice shading from the selected color's hue
    h = LCHuv(parse(Colorant, basecolor)).h
    c = sequential_palette(h, 256)

    #Read the simulation data for each size
    i=1
    for n in system_sizes
        #read_data!(data_path, n, nparts, ngammas, nsims, av_r, av_sus, gammas)
        read_data_simpler!("$data_path/diagrams_$n", ngammas, nsims, av_r, av_sus, gammas)

        #Plot the simulation data
        scatter!(axis, gammas, av_r, color=fscolors[i], label=L"N=%$n")

        #Plot the theory line
        if use_OA
            r_th = TheoryFormulas.finite_size_r(gammas, s2, n)
            lines!(axis, gammas, r_th, color=fscolors[i])
        else
            q_th, r_th = TheoryFormulas.diagram_amplitudes(6, n, 0.0, 0.2, ngammas, s2, ""; writetofile=false)
            lines!(axis, q_th, r_th, color=fscolors[i])
        end
        i += 1
    end

    #Axis label and legend
    axis.xlabel = L"J"
    axislegend(axis, position=(0.005, 0.95))
end

#Start the figure with two axes 
#set_theme!(StyleFuncs.two_col_figure(3.6))
set_theme!(StyleFuncs.one_col_figure(1.5))
fig = Figure(figure_padding=5)
group = fig[1,1] = GridLayout()
axs = [Axis(group[j,1], xgridvisible=false, ygridvisible=false) for j=1:2]

#Gap between plots 
#colgap!(group, 10)
#rowgap!(group, 10)

#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [colors[i] for i in [2,3,1]]

#Plot in each axes
plot_thermodynamic_limit(axs[1], colors, 0.1, "kuramoto_julia_bien"; ngammas=43, eqnumber="Eq. (9)")
plot_finite_sizes(axs[2], colors[3])

#Fine tuning of the decorations
#hideydecorations!(axs[2])
#for ax in axs
#    ax.xticks = [0.0, 0.1, 0.2] 
#    #ylims!(ax, 0.0, 0.85)
#end

hidexdecorations!(axs[1])
for ax in axs
    ax.yticks = [0.0, 0.5, 0.8] 
#    #ylims!(ax, 0.0, 0.85)
end

#Add axes labels
StyleFuncs.label_axes(axs, pos=[0.05, 0.85])

#Save the figure
save("figure1_alt.pdf", fig, pt_per_unit = 1)

# ------

#Supplementary figure
set_theme!(StyleFuncs.two_col_figure(3.5))
fig = Figure(figure_padding=7)

axs = [Axis(fig[1,1], title=L"N=10^5, \quad \sigma^2 = 0.5"), Axis(fig[1,2], title=L"Finite-size Ott-Antonsen eq. (D10), $\sigma^2 = 0.1$") ]

plot_thermodynamic_limit(axs[1], colors, 0.5, "kuramoto_julia_bien_s05"; ngammas=43, nsims=30, eqnumber="Main text eq. (9)")
plot_finite_sizes(axs[2], colors[1]; use_OA = true)

StyleFuncs.label_axes(axs)

axs[2].xticks = [0.0, 0.1, 0.2] 
save("figure1_sup.pdf", fig, pt_per_unit = 1)


# -------------

# Larger format, separate figs, adequate for talks
set_theme!(StyleFuncs.slide_figure(2.))
fig = Figure(figure_padding=5)

#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [colors[i] for i in [2,3,1]]

#Plot in each axes
ax = Axis(fig[1,1])
plot_thermodynamic_limit(ax, colors, 0.1, "kuramoto_julia_bien"; ngammas=43, eqnumber="Eq. (9)")
ax.yticks = [0.0, 0.5, 0.8] 


#Save the figure
save("figure_talks_phasediagrams.pdf", fig, pt_per_unit = 1)

fig = Figure(figure_padding=5)

#Plot in each axes
ax = Axis(fig[1,1])
plot_finite_sizes(ax, colors[3])
ax.yticks = [0.0, 0.5, 0.8] 
ax.ylabel = L"R"

#Save the figure
save("figure_talks_fs.pdf", fig, pt_per_unit = 1)