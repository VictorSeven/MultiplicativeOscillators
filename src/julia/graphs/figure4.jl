
using CairoMakie
using DelimitedFiles
using LsqFit

include("style_funcs.jl")

function plot_variance(axis)
    nparts = 3
    ngammas = 50
    nsims = 100
    n = 100000
    s2 = 0.1
    data_path = "../../../data/diagrams"

    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas)


    for index=1:nsims
        data = Matrix{Float64}(undef, 0, 4) 
        for part=1:nparts
            nextdata = readdlm("$data_path/diagrams_$n/diagram_$(index-1)_$(part-1)")
            data = vcat(data, nextdata) 
        end
        av_r += data[:, 2] 
        av_sus += data[:, 3] 
        gammas = data[:, 1]
    end
    av_r /= nsims
    av_sus /= nsims

    var_theo = readdlm("../../../data/diagrams/theoretical/stochastic_full_50")[:,2] #integrate_full.(0.1, gammas, s2)
    teogammas = LinRange(0.0, 0.2, 100)
    #scatter!(axis, teogammas, 2*var_theo, label="Full System", markersize=3)
    lines!(axis, teogammas, var_theo, label="Full System")
    scatter!(axis, gammas, av_sus, markersize=4, color=:gray, label="Simulation")


    create_legend(axis, (-6.5, 2))

    axis.xlabel = L"J"
    axis.ylabel = L"\chi"

    xlims!(axis, 0.095, 0.125)

end

function plot_gamma(axis)
    nparts = 3
    ngammas = 50
    nsims = 100
    n = 100000
    s2 = 0.1
    data_path = "../../../data/diagrams"

    av_r = zeros(ngammas) 
    av_sus = zeros(ngammas) 
    gammas = Vector{Float64}(undef, ngammas)


    for index=1:nsims
        data = Matrix{Float64}(undef, 0, 4) 
        for part=1:nparts
            nextdata = readdlm("$data_path/diagrams_$n/diagram_$(index-1)_$(part-1)")
            data = vcat(data, nextdata) 
        end
        av_r += data[:, 2] 
        av_sus += data[:, 3] 
        gammas = data[:, 1]
    end
    av_r /= nsims
    av_sus /= nsims

    var_theo = readdlm("../../../data/diagrams/theoretical/stochastic_full_50")[:,2] #integrate_full.(0.1, gammas, s2)
    teogammas = LinRange(0.0, 0.2, 100)
    #var_theo = av_sus

    gc = 0.11

    epsilon = (teogammas .- gc) / gc
    #epsilon = (gammas .- gc) / gc


    logvar = log.(abs.(var_theo[epsilon .> 0]))
    #var_theo = log.(abs.(av_sus[epsilon .> 0]))
    logepsilon = log.(epsilon[epsilon .> 0])
    println(epsilon[epsilon .> 0])
    println(epsilon[epsilon .< 0])

    p=10
    logvar= logvar[begin:begin+p]
    logepsilon = logepsilon[begin:begin+p] 


    @. model(x,p) = p[1]*x+p[2] 
    fit = curve_fit(model, logepsilon, logvar, [1.0, -10])
    b = fit.param
    sigma = margin_error(fit)
    println(b, " ", sigma)

    lines!(axis[1], logepsilon, model(logepsilon, b), color=:black)
    scatter!(axis[1], logepsilon, logvar, label="Full System", markersize=3)


    gc = 0.102

    epsilon = (teogammas .- gc) / gc
    #epsilon = (gammas .- gc) / gc

    #logvar= log.(abs.(av_sus[epsilon .< 0]))
    logvar = log.(abs.(var_theo[epsilon .< 0]))
    logepsilon = log.(-epsilon[epsilon .< 0])

    logvar = logvar[end-p:end]
    logepsilon = logepsilon[end-p:end] 

    @. model(x,p) = p[1]*x+p[2] 
    fit = curve_fit(model, logepsilon, logvar, [1.0, -10])
    b = fit.param
    sigma = margin_error(fit)
    println(b, " ", sigma)

    lines!(axis[2], logepsilon, model(logepsilon, b), color=:black)
    scatter!(axis[2], logepsilon, logvar, label="Full System", markersize=3)


    #=
    scatter!(axis, epsilon, log.(var_theo), label="Full System", markersize=3)
    #ylims!(axis, 0, 1.5)
    #scatter!(axis, gammas, av_sus, markersize=4, color=:gray, label="Simulation")

    epsilon = LinRange(0.0, 1.0, 10)
    teoline =  epsilon .- 12 
    lines!(axis, epsilon , teoline, color=:black)


    #create_legend(axis, (-6.5, 2))

    axis.xlabel = L"J"
    axis.ylabel = L"\chi"
    =#

end

fig = Figure(resolution=one_col_size(2.25), fontsize=9, figure_padding=7)
group = fig[1,1] = GridLayout()
left  = group[1,1] = GridLayout()
right = group[1,2] = GridLayout()
axvar = Axis(left[1,1], xgridvisible=false, ygridvisible=false) 
axgamma = [Axis(right[j,1], xgridvisible=false, ygridvisible=false) for j=1:2]

axs = [axvar, axgamma[1], axgamma[2]]

colgap!(group,5)
rowgap!(right, 2)


hidespines!.(axs, :t, :r)

linkxaxes!(axgamma[1], axgamma[2])

hideydecorations!.(axs, label=true)
hidexdecorations!(axgamma[1])

plot_variance(axvar)
plot_gamma(axgamma)

label_axes(axs)
save("figure4.pdf", fig, pt_per_unit = 1)