using CairoMakie
using DelimitedFiles
using LsqFit
using Statistics

include("style_funcs.jl")
include("palette.jl")

using .ArtsyPalettes
using .StyleFuncs


function read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, gammas; ncols=4)
    #Initialize vectors
    av_r   .= zeros(ngammas) 
    av_sus .= zeros(ngammas) 
    gammas .= Vector{Float64}(undef, ngammas)

    #Read simulation data. There are several simulations
    for index=0:nsims-1
        data = readdlm("$data_path/diagram_sim$(index)")

        #Get ensemble averages over all simulations 
        av_sus .+= data[1:ngammas, 3] 
        gammas .= data[1:ngammas, 1]
        
    end

    #Finish means
    av_r ./= nsims
    av_sus ./= nsims
end

#TODO resolver esto
function numerical_derivative(y, x)
    derivs = Vector{Float64}(undef, length(x))

    h_next = x[begin+2:end-1] - x[begin+1:end-2]
    h_back = x[begin:end] - x[begin-1:end-1]

    f_next = y

    @. derivs[begin+1:end-1] = (y[begin+2:end] - y[begin:end-2]) / (x[begin+2:end] - x[begin:end-2])  

    derivs[1] = (y[2] - y[1]) / (x[2] - x[1])  
    derivs[end] = (y[end] - y[end-1]) / (x[end] - x[end-1])  

    return derivs


    #return @views @. (y[2:end] - y[1:end-1]) / (x[2:end] - x[1:end-1]) 
end

function plot_effective_exponent(ax, q, susc)
    peak, halfpoint = findmax(susc) 
    rc = q[halfpoint] 
    eps = @. (q - rc) / rc 
    gammaeff = numerical_derivative(log.(susc), log.(abs.(eps)))

    lin2fit(x, p) = p[1]*x .+ p[2] 
    p0 = [1., 0.]
    offset = 1


    x = eps[begin:halfpoint-1]
    y = -1 ./ gammaeff[begin:halfpoint-1]
    scatter!(ax, x, y)

    println(x[begin:end-offset])
    fit = curve_fit(lin2fit, x[begin:end-offset], y[begin:end-offset], p0)
    p, err = fit.param, estimate_errors(fit, 0.95)
    println("$(p[1]) ± $(err[1])")
    println("$(p[2]) ± $(err[2])")
    lines!(ax, x, lin2fit(x, p), label="γ = $(round(p[2], sigdigits=3)) ± $(round(err[2], sigdigits=1))")

    x = eps[halfpoint+1:end]
    y = -1 ./ gammaeff[halfpoint+1:end]
    scatter!(ax, x, y)

    fit = curve_fit(lin2fit, x[begin+offset:end], y[begin+offset:end], p0)
    p, err = fit.param, estimate_errors(fit, 0.95)
    println("$(p[1]) ± $(err[1])")
    println("$(p[2]) ± $(err[2])")
    lines!(ax, x, lin2fit(x, p), label="γ' = $(round(p[2], sigdigits=3)) ± $(round(err[2], sigdigits=1))")


    axislegend(ax, position=(0.1, 0.8))
    vlines!(ax, [0.0], color=:black)


    
    xlims!(ax, -0.45, 0.45)
    ylims!(ax, 0.5, 1.5)

    ax.xlabel =  L"$\varepsilon$"
    ax.ylabel =  L"$\gamma_{\text{eff}}$"

end


#Start the figure with two axes 
#fig = Figure(resolution=two_col_size(2*1.618), fontsize=9, figure_padding=7)
set_theme!(two_col_figure(2*1.618))
fig = Figure(figure_padding=7)

ax = Axis(fig[1,1])

ax = Axis(fig[1,2])

ngammas = 100
nsims = 3000 
data_path = "../../../data/gamma/offcritical_n1e6_sigsq0.1"
s2 = 0.1

#Intialize vectors
av_r = Vector{Float64}(undef, ngammas)
av_sus = Vector{Float64}(undef, ngammas)
gammas = Vector{Float64}(undef, ngammas)

read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, gammas)

#Get the colormap
colors = met_brew("Egypt")
colors = [colors[i] for i in [2,3]]

plot_effective_exponent(ax, gammas, av_sus)



#Save the figure
save("figure2.pdf", fig, pt_per_unit = 1)