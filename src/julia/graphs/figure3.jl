using CairoMakie
using Distributions

include("style_funcs.jl")

function get_corrs(x, y, a=0.01)
    corr = zeros(2,2)
    corr[1,1] = a * (1.0 - x^2 + y^2)
    corr[1,2] = -2*a*x*y
    corr[2,1] = corr[1,2]
    corr[2,2] = a * (1.0 + x^2 - y^2)
    return corr
end

fig = Figure(resolution=one_col_size(1), fontsize=9, figure_padding=7)
axs = Axis(fig[1,1], xgridvisible=false, ygridvisible=false) 

hidespines!.(axs, :t, :r, :b, :l)

hlines!(0, color=:black)
vlines!(0, color=:black)

amplitudes = [0, 0.2, 0.6, 0.95]
angles = LinRange(0.0, 7π/4, 8)


coords = Matrix{Float64}(undef, 4*8, 2)
index = 1

blacks = to_colormap(:grayC) 
blacks_used = [blacks[i] for i in [1, 60, 100, 150]]

for (i,r) in enumerate(amplitudes)
    arc!(axs, Point2f(0.0, 0.0), r, 0, 2π, color=blacks_used[i])
    
    for phi in angles
        x = r*cos(phi)
        y = r*sin(phi)
        corr = get_corrs(x,y)

        mvn = MvNormal([x,y], corr) 
        xi = rand(mvn, 100)

        rand_angles = atan.(xi[2,:], xi[1,:])
        rand_angles = mod2pi.(rand_angles) 

        scatter!(axs, xi[1,:], xi[2, :], markersize=2, color=rand_angles, colormap=:vikO, colorrange=(0,2π))

        coords[index, :] = [x,y]
        global index += 1
    end
end


epslim = 0.05
xlims!(axs, -1-epslim, 1+epslim)
ylims!(axs, -1-epslim, 1+epslim)

arc!(axs, Point2f(0.0, 0.0), 1, 0, 2π, color=:black)
text!(axs, 0.47, 0.85; text=L"|Z_1|=1")

hidexdecorations!(axs)
hideydecorations!(axs)

save("figure3.pdf", fig, pt_per_unit = 1)