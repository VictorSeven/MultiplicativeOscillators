using DelimitedFiles
using Statistics
using CairoMakie

include("style_funcs.jl")

function angle_diff(angle1, angle2)
    angle1_mod = mod.(angle1, 2π)
    angle2_mod = mod.(angle2, 2π)

    diff = @. abs(angle1_mod - angle2_mod) 
    return @. min(diff, 2π - diff) 
end

function diff_random_angles(n_angles, its, nharms)
    k=collect(1:nharms)
    av_diff = zeros(nharms)
    for i=1:its 
        phis = 2π*rand(n_angles)
        zk = [mean(exp.(im*phis*k)) for k=1:nharms]
        psi = angle.(zk)
        av_diff += angle_diff.(psi, k*psi[1])
    end
    av_diff /= its
    return mod.(av_diff, 2π)
end

function diff_gaussian_angles(n_angles, its, nharms)
    k=collect(1:nharms)
    av_diff = zeros(nharms)
    for i=1:its 
        phis = 2*π*rand() .+ 0.5 * randn(n_angles)
        zk = [mean(exp.(im*phis*k)) for k=1:nharms]
        psi = angle.(zk)
        av_diff += angle_diff.(psi, k*psi[1])
    end
    av_diff /= its
    return mod.(av_diff, 2π)
end

function analyse_angle_recurrence(data_path, ax, label)

    xy = readdlm(data_path)
    nharms = (size(xy)[2] - 1) ÷ 2

    z    = xy[:,2] .+ 1.0im .* xy[:,3]
    psi1 = angle.(z)

    psidiff = Vector{Float64}(undef, nharms)
    psidiff[1] = 0.0

    for k=2:nharms
        z = xy[:,2*k] .+ 1.0im .* xy[:,2*k+1]
        psik = angle.(z)
        psidiff[k] = mean(angle_diff(psik, k*psi1))
    end

    k = collect(1:nharms) 
    scatter!(ax, k, psidiff, markersize=5, label=label)

    ylims!(ax, 0.0, π)
end

fig = Figure(resolution=one_col_size(2), fontsize=9, figure_padding=5)
ax = Axis(fig[1,1], xlabel=L"k", ylabel=L"\left| \psi_k - k\psi_1 \right|", 
xgridvisible=false, ygridvisible=false)

hidespines!(ax, :t, :r)

#data_path = "../../../data/timeseries/series_micro/series_100/series_0.13_1" 
#data_path = "../../../data/timeseries/time_micro_asyn_100_4"
data_path = "../../../data/series_micro/size_100/series_long_0.07"
analyse_angle_recurrence(data_path, ax, "Kuramoto")

#data_path = "../../../data/timeseries/time_micro_asyn_EXC_100_1"
#analyse_angle_recurrence(data_path, ax, "Excitable Osc.")

#scatter!(ax, k, diff_random_angles(100, 10000, 10), markersize=5, label="Random")
#scatter!(ax, k, diff_gaussian_angles(10000, 10000, 10), markersize=5, label="Gaussian")

data_path = "../../../data/series_micro/size_100/series_long_0.2"
analyse_angle_recurrence(data_path, ax, "Synchro")

data_path = "../../../data/series_micro/size_10000/series_long_0.2"
analyse_angle_recurrence(data_path, ax, "Synchro")

hlines!(ax, π/2, color=:gray, linestyle=:dash, label=L"\pi/2")

create_legend(ax, (0.005, 1.2))

save("figure2.pdf", fig, pt_per_unit=1)
