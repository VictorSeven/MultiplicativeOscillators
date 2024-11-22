using CairoMakie
using Distributions
using LinearAlgebra

include("style_funcs.jl")

using .StyleFuncs

function get_corrs(x, y, a=0.005)
    corr = zeros(2,2)
    corr[1,1] = a * (1.0 - x^2 + y^2)
    corr[1,2] = -2*a*x*y
    corr[2,1] = corr[1,2]
    corr[2,2] = a * (1.0 + x^2 - y^2)
    return corr
end

set_theme!(StyleFuncs.one_col_figure(1.15)) #perfect circle: 1.05
fig = Figure(fontsize=9, figure_padding=3)
grid = fig[1,1] = GridLayout()
axs = Axis(grid[1,1])
insets = [Axis(grid[2,1]), Axis(grid[1,2], rightspinevisible=true)]
insets[2].yaxisposition = :right

rowsize!(grid, 1, Relative(0.75))
colsize!(grid, 1, Relative(0.75))
rowgap!(grid, 1.)
colgap!(grid, 1.)

hidespines!.(axs, :t, :r, :b, :l)
hidexdecorations!(axs)
hideydecorations!(axs)

hidespines!.(insets[1], :t, :r, :l)
#hideydecorations!(insets[1])
insets[1].yticklabelsvisible = false
insets[1].yticksvisible = false
hidespines!.(insets[2], :t, :b, :l)
hidexdecorations!(insets[2])

amplitudes = [0.4, 0.7, 0.975]
angles = LinRange(0.0, 7π/4, 8)

colormap = :hawaii
cmap =  cgrad(colormap)
clrs = [cmap.colors[i] for i=1:256÷4:256]

corr = get_corrs(0., 0.)
mvn = MvNormal([0., 0.], corr) 
xi = rand(mvn, 500)
scatter!(axs, xi[1,:], xi[2, :], markersize=2, color=clrs[1])

for (i,r) in enumerate(amplitudes)
    arc!(axs, Point2f(0.0, 0.0), r, 0, 2π, color=clrs[1+i], linewidth=0.5)

    for (j,phi) in enumerate(angles)
        x = r*cos(phi)
        y = r*sin(phi)
        corr = get_corrs(x,y)

        mvn = MvNormal([x,y], corr) 
        xi = rand(mvn, 300)

        scatter!(axs, xi[1,:], xi[2, :], markersize=2, color=clrs[1 + i])
    end
end



epslim = 0.05
xlims!(axs, -1-epslim, 1+epslim)
ylims!(axs, -1-epslim, 1+epslim)

arc!(axs, Point2f(0.0, 0.0), 1, 0, 2π, color=:black)
text!(axs, 0.47, 0.85; text=L"R_1=1", fontsize=12)



x = LinRange(-1, 1, 100)
fluc = @. sqrt(1 -x^2)

clrs = cmap[abs.(x)]
println(clrs)
lines!(insets[1], x, fluc, color=clrs, linewidth=2)
insets[1].xlabel = L"x"
insets[1].ylabel = L"\langle \xi_{R_1} ^2  \rangle"
insets[1].xticks = [-1, 0, 1]

lines!(insets[2], -fluc, x, color=clrs, linewidth=2)
insets[2].ylabel = L"y"
insets[2].yticks = [-1, 0, 1]

save("figure3.pdf", fig, pt_per_unit = 1)