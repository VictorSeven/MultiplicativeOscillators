
using CairoMakie
using DelimitedFiles
using KernelDensity
using Printf

include("style_funcs.jl")
include("palette.jl")

using .ArtsyPalettes
using .StyleFuncs


"""
Plot the distribution of Zk obtained with coupling q in the selected axis.
The data is inside data_path and the user can select the bounds of the bins and colors
"""
function plot_distributions!(q, axis, data_path, bounds, colors)

    filenames = ["additive_noFS", "additive", "amplitude", "kuramoto"]
    labels = ["Additive", "Add + FS drift", "Eqs. (8)", "Simulation"]

    q2d = @sprintf("%.2f", q)

    binspace = bounds[1]:(bounds[2]-bounds[1])/40:bounds[2]

    for i in eachindex(filenames)
        name = filenames[i]
        c = colors[i]
        lab = labels[i]
        rsamples = readdlm("$(data_path)/series_$(name)_$(q2d)")
        stephist!(axis, abs.(rsamples[:,2]), bins=binspace, color=c, label=lab, normalization=:pdf)
    end


    axis.xticks = [bounds[1], bounds[2]]
    #axis.yticks = [axis.yticks[1], axis.yticks[end]]
    axis.title = "J=$q"
    axis.titlefont = :regular
end


set_theme!(StyleFuncs.one_col_figure(3.0))
fig = Figure(figure_padding=(2, 7, 2, 1), backgroundcolor=:transparent)
group = fig[1,1] = GridLayout()
axs = [Axis(group[2,j], xgridvisible=false, ygridvisible=false) for j=1:3]

hidespines!.(axs, :t, :r)

data_path = "../../../data/series4dists/series_1000"

colors = ArtsyPalettes.met_brew("Isfahan1")
colors = [colors[i] for i in [2, 4, 6]]
append!(colors, ["#000000"])
plot_distributions!(0.05, axs[1], data_path, (0., 0.2), colors)
plot_distributions!(0.1, axs[2], data_path, (0., 0.5), colors)
plot_distributions!(0.2, axs[3], data_path, (0.75, 0.9), colors)

axs[2].xlabel = "R"
axs[2].xlabelpadding = -4.
axs[1].ylabel = "p(R)"

leg = Legend(group[1,1:3], axs[1], position=(0., 0.), orientation=:horizontal, colgap=5, backgroundcolor=:transparent)
#Label(group[2, 1, Top()], "Title")


rowgap!(group, 0)
colgap!(group, 3)
rowsize!(group, 2, Relative(0.7))

StyleFuncs.label_axes(axs, pos=[0.75, 0.8, 0.8])

save("figure4.pdf", fig, pt_per_unit = 1)



# ----- Same figure, larger format

set_theme!(StyleFuncs.slide_figure(2.))
fig = Figure(backgroundcolor=:transparent)
group = fig[1,1] = GridLayout()
axs = [Axis(group[2,j], xgridvisible=false, ygridvisible=false) for j=1:3]

hidespines!.(axs, :t, :r)

data_path = "../../../data/series4dists/series_1000"

colors = ArtsyPalettes.met_brew("Isfahan1")
colors = [colors[i] for i in [2, 4, 6]]
append!(colors, ["#000000"])
plot_distributions!(0.05, axs[1], data_path, (0., 0.2), colors)
plot_distributions!(0.1, axs[2], data_path, (0., 0.5), colors)
plot_distributions!(0.2, axs[3], data_path, (0.75, 0.9), colors)

axs[2].xlabel = "R"
axs[2].xlabelpadding = -4.
axs[1].ylabel = "p(R)"

leg = Legend(group[1,1:3], axs[1], position=(0., 0.), orientation=:horizontal, colgap=5, backgroundcolor=:transparent)


rowgap!(group, 0)
colgap!(group, 3)
rowsize!(group, 2, Relative(0.7))


save("figure4_talks.pdf", fig, pt_per_unit = 1)
