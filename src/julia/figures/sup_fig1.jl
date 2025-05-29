using CairoMakie
using DelimitedFiles

include("style_funcs.jl")
include("palette.jl")

using .StyleFuncs
using .ArtsyPalettes

"""
Distance between 2 angles in [0, 2π]
"""
function angledist(a1, a2)
    return min(2π - abs(a1 - a2), abs(a1- a2))
end

"""
Add the nice grey color background to the figure
"""
function background_fig(ax, start_harm, nharm; sqham=4, c1="#dfdfdf", c2="bfbfbf")
    vspan!(ax, start_harm, nharm, color="#dfdfdf", strokewidth=0.)
end

"""
Plot the figure in two selected axs, by taking the data from filename, 
using c1 and c2 as starting colors, and start_harm as the first harmonic that will be considered not valid 
"""
function plot_figure(axs, filename, c1, c2, start_harm)

    #Constants
    nreps = 100
    nharm = 20 
    N = 100000

    # --- All simulations 
    #Do the logscale by hand so the errors are shown between the logscale curves

    ravg = zeros(nharm)
    logravg = zeros(nharm)
    logrstd = zeros(nharm)

    davg = zeros(nharm-1)
    dstd = zeros(nharm-1)

    #Average over all simulations
    for i=0:nreps-1
        zk = readdlm("$(filename)_$i")
        zkcomplex = zk[1:2:2*nharm] .+ 1im * zk[2:2:2*nharm]
        r = abs.(zkcomplex)
        psi = angle.(zkcomplex)
        logr = log10.(r)
        ravg .+= r

        logravg .+= logr
        logrstd .+= logr .^2

        angled = angledist.(psi[1:end-1], psi[2:end]) .- abs(psi[1]) 
        davg .+= angled 
        dstd .+= angled .^2
    end

    #Finish averages
    ravg ./= nreps
    logravg ./= nreps
    davg ./= nreps
    @. logrstd = logrstd / nreps - logravg^2
    @. dstd = dstd / nreps - davg^2


    #Plot
    ax = axs[1] 
    background_fig(ax, start_harm, nharm)
    hlines!(ax, log10(1/N), color=:gray, linestyle=:dash)

    plot!(ax, 1:nharm, logravg, color=c1)
    lines!(ax, 1:nharm, logravg, color=c1, label="Simulation")
    errorbars!(ax, 1:nharm, logravg, logrstd, color=c1)
    lines!(ax, 1:nharm, (1:nharm) .* log10.(ravg[1]), color=c2, label="Ott-Antonsen")
    plot!(ax, 1:nharm, (1:nharm) .* log10.(ravg[1]), color=c2)


    ax = axs[2]
    background_fig(ax, start_harm, nharm)
    hlines!(ax, 0., color=c2)

    plot!(ax, 1:nharm-1, davg, color=c1)
    errorbars!(ax, 1:nharm-1, davg, dstd, color=c1)
    ylims!(ax, -pi/2, pi/2)
end


#Start the figure with four axes 
set_theme!(StyleFuncs.two_col_figure(1.5))
fig = Figure(figure_padding=5)
group = fig[1,1] = GridLayout()

axs = [Axis(group[1,1])]
push!(axs, Axis(group[1,2]))
push!(axs, Axis(group[2,1]))
push!(axs, Axis(group[2,2]))


#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [:black, colors[2]]

#Set titles and labels
axs[1].title = "J = 0.12"
axs[2].title = "J = 0.3"

axs[1].ylabel = L"$R_k"
axs[3].ylabel = L"$\psi _{k+1} - (\psi_k + \psi_1)$"

#Set the correct ticks
for axi in [1,2]
    ax = axs[axi]
    ylims!(ax, -6, 1)
    ax.yticks = ([-k for k=0:7], [rich("10", superscript("$k")) for k=0:7])
end

for axi in [3,4]
    ax = axs[axi]
    ax.yticks = ([-pi/2, -pi/4, 0, pi/4, pi/2], ["-π/2", "-π/4", "0", "π/4", "π/2"])
    ax.xlabel = "Harmonic k"
end

#Do the plot for each data file
plot_figure([axs[1], axs[3]], "../../../data/ansatz/30harms_q012",  colors[1], colors[2], 5)
plot_figure([axs[2], axs[4]], "../../../data/ansatz/30harms_q03", colors[1], colors[2], 10)

axislegend(axs[1], position=(0.95, 0.8))

#Add axes labels
StyleFuncs.label_axes(axs, pos=[0.05, 0.9])

#Save the figure
set_theme!(StyleFuncs.two_col_figure(3.5))
save("figure2_sup.pdf", fig, pt_per_unit = 1)
