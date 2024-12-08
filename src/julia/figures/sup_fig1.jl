using CairoMakie
using DelimitedFiles

include("style_funcs.jl")
include("palette.jl")

using .StyleFuncs
using .ArtsyPalettes

function angledist(a1, a2)
    return min(2π - abs(a1 - a2), abs(a1- a2))
end

function background_fig(ax, start_harm, nharm; c1="#dfdfdf", c2="bfbfbf")
    #vspan!(ax, from, to, color="#dfdfdf", strokewidth=0.)
    vspan!(ax, start_harm, nharm, color="#cfcfcf", strokewidth=0.)
end

function plot_figure(axs, filename, c1, c2)
    nreps = 100
    nharm = 20 
    start_harm = 11
    N = 100000

    # --- Single simulation

    zk = readdlm("$(filename)_0")
    zkcomplex = zk[1:2:2*nharm] .+ 1im * zk[2:2:2*nharm]
    r = abs.(zkcomplex)
    psi = angle.(zkcomplex)

    ax = axs[1] 
    background_fig(ax, start_harm, nharm)
    #hlines!(ax, [1/N, 1/sqrt(N)], color=:gray, linestyle=:dash)
    hlines!(ax, 1/N, color=:gray, linestyle=:dash)

    lines!(ax, 1:nharm, r, color=c1, label="Simulation")
    plot!(ax,  1:nharm,  r, color=c1)
    lines!(ax,  1:nharm,  r[1] .^ (1:nharm), color=c2, label="Ott-Antonsen")
    plot!(ax,  1:nharm,  r[1] .^ (1:nharm), color=c2)

    ax = axs[2] 
    background_fig(ax, start_harm, nharm)
    hlines!(ax, 0., color=c2)

    plot!(ax, 1:nharm-1, angledist.(psi[1:end-1], psi[2:end]) .- abs(psi[1]), color=c1)
    ylims!(ax, -pi/2, pi/2)


    # --- All simulations 
    #Do the logscale by hand so the errors are shown between the logscale curves

    ravg = zeros(nharm)
    logravg = zeros(nharm)
    logrstd = zeros(nharm)

    davg = zeros(nharm-1)
    dstd = zeros(nharm-1)

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

    ravg ./= nreps
    logravg ./= nreps
    davg ./= nreps
    @. logrstd = logrstd / nreps - logravg^2
    @. dstd = dstd / nreps - davg^2


    fig = Figure()
    ax = axs[3] 
    background_fig(ax, start_harm, nharm)
    #hlines!(ax, log10.([1/N, 1/sqrt(N)]), color=:gray, linestyle=:dash)
    hlines!(ax, log10(1/N), color=:gray, linestyle=:dash)

    plot!(ax, 1:nharm, logravg, color=c1)
    lines!(ax, 1:nharm, logravg, color=c1)
    errorbars!(ax, 1:nharm, logravg, logrstd, color=c1)
    lines!(ax, 1:nharm, (1:nharm) .* log10.(r[1]), color=c2)
    plot!(ax, 1:nharm, (1:nharm) .* log10.(r[1]), color=c2)

    ax = axs[4]
    background_fig(ax, start_harm, nharm)
    hlines!(ax, 0., color=c2)

    plot!(ax, 1:nharm-1, davg, color=c1)
    errorbars!(ax, 1:nharm-1, davg, dstd, color=c1)
    ylims!(ax, -pi/2, pi/2)


end


#Start the figure with two axes 
#set_theme!(StyleFuncs.two_col_figure(3.6))
set_theme!(StyleFuncs.two_col_figure(1.5))
fig = Figure(figure_padding=5)
group = fig[1,1] = GridLayout()
axs = [Axis(group[1,1], yscale=log10)]
push!(axs, Axis(group[1,2]))
push!(axs, Axis(group[2,1]))
push!(axs, Axis(group[2,2]))
#axs = [Axis(group[j,1], xgridvisible=false, ygridvisible=false) for j=1:2]

#Gap between plots 
#colgap!(group, 10)
#rowgap!(group, 10)

#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [:black, colors[1]]


for axi in [1,3]
    ax = axs[axi]
    ax.xlabel = "Harmonic k"
    ax.ylabel = L"$R_k"
end
for axi in [2,4]
    ax = axs[axi]
    ax.ylabel = L"$\psi _{k+1} - (\psi_k + \psi_1)$"
    ax.yticks = ([-pi/2, -pi/4, 0, pi/4, pi/2], ["-π/2", "-π/4", "0", "π/4", "π/2"])
end
for axi in [3,4]
    ax = axs[axi]
    ax.xlabel = "Harmonic k"
end

plot_figure(axs, "../../../data/ansatz/30harms_q012", colors[1], colors[2])
axislegend(axs[1], position=(0.95, 0.8))

#Add axes labels
StyleFuncs.label_axes(axs, pos=[0.05, 0.85])

#Save the figure
set_theme!(StyleFuncs.two_col_figure(3.5))
save("figure1_sup.pdf", fig, pt_per_unit = 1)
