using CairoMakie
using DelimitedFiles
using LsqFit
using Statistics

include("style_funcs.jl")
include("palette.jl")

using .ArtsyPalettes
using .StyleFuncs


"""
Reads the data at data_path to obtain mean and variance stored in the files
"""
function read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, var_r, gammas)

    #Initialize vectors
    av_r   .= 0. 
    var_r  .= 0. 
    av_sus .= 0. 

    #Read simulation data. There are several simulations
    for index=0:nsims-1
        data = readdlm("$data_path/diagram_sim$(index)")

        #Get ensemble averages over all simulations 
        av_sus .+= data[1:ngammas, 3] 
        var_r  .+= data[1:ngammas, 2] .* data[1:ngammas, 2] 
        av_r   .+= data[1:ngammas, 2] 

        gammas .= data[1:ngammas, 1]
        
    end

    #Finish means
    av_r ./= nsims
    var_r./= nsims
    @. var_r = var_r - av_r * av_r
    av_sus ./= nsims
end

"""
Compute the derivative of the function y(x) which has been sampled in discrete steps. The implementation here has
been copied from Numpy's gradient function (centered finite differences) with boundary conditions having first order.
"""
function numerical_derivative(y, x)

    #Create target vector
    n = length(x)
    derivs = Vector{Float64}(undef, n)

    #Take care of the boundary conditions
    derivs[1] = (y[2] - y[1]) / (x[2] - x[1]) 
    derivs[end] = (y[end] - y[end-1]) / (x[end] - x[end-1]) 

    #dx for non-equally spaced points
    xdiff = @views x[2:end] - x[1:end-1]
    dx1 = xdiff[1:end-1]
    dx2 = xdiff[2:end]

    #Constants for left, centered and right steps
    a = @. -(dx2)/(dx1 * (dx1 + dx2))
    b = @. (dx2 - dx1) / (dx1 * dx2)
    c = @. dx1 / (dx2 * (dx1 + dx2))

    #Slices 
    slice1 = 2:n-1
    slice2 = 1:n-2
    slice3 = 2:n-1
    slice4 = 3:n
    
    #Finish computation and return values
    @. derivs[slice1] = @views a*y[slice2] + b*y[slice3] + c*y[slice4]  

    return derivs
end

"""
Plot the effective exponent at axis ax via the definition of critical exponent, dlog(susc)/dlogq 
"""
function plot_effective_exponent(ax, q, susc)

    #Get the maximum of the susceptibility
    peak, halfpoint = findmax(susc) 
    rc = q[halfpoint] 
    eps = @. (q - rc) / rc 
    #Compute effective value of the exponent
    gammaeff = numerical_derivative(log.(susc), log.(abs.(eps)))
    offset = 1 #To avoid taking the nan values at criticality

    #Prepare functions for fitting
    lin2fit(x, p) = p[1]*x .+ p[2] 
    quad2fit(x, p) = p[1]*x .+ p[3]*x .^2 .+ p[2] 

    #Do the scatterplot
    x = eps[begin:halfpoint-offset]
    y = -1 ./ gammaeff[begin:halfpoint-offset]
    scatter!(ax, x, y)

    #Do fit for gamma and plot it
    p0 = [0., 1.]
    fit = curve_fit(lin2fit, x[begin:end-offset], y[begin:end-offset], p0)
    p, err = fit.param, estimate_errors(fit, 0.68)
    #println("$(p[1]) ± $(err[1])")
    #println("$(p[2]) ± $(err[2])")
    lines!(ax, x, lin2fit(x, p), label="γ = $(round(p[2], sigdigits=4)) ± $(round(err[2], sigdigits=1))")


    #Scatterplot for the remaining things
    x = eps[halfpoint+1:end]
    y = -1 ./ gammaeff[halfpoint+1:end]
    scatter!(ax, x, y)

    #Fit for gamma'...
    p0 = [0., 1., 0.]
    fit = curve_fit(quad2fit, x[begin+offset:end], y[begin+offset:end], p0)
    p, err = fit.param, estimate_errors(fit, 0.68)

    lines!(ax, x, quad2fit(x, p), label="γ' = $(round(p[2], sigdigits=4)) ± $(round(err[2], sigdigits=1))")

    #Make the axes look nice
    vlines!(ax, [0.0], color=:black)
    
    xlims!(ax, -0.45, 0.45)
    ylims!(ax, 0.5, 1.5)

    ax.xlabel =  L"$\varepsilon$"
    ax.ylabel =  L"$\gamma_{\text{eff}}$"

end


#Start the figure with two axes 
set_theme!(StyleFuncs.one_col_figure(1.5))
fig = Figure(figure_padding=2)
axs = [Axis(fig[j,1]) for j=1:2]

#Initialize the vectors we'll need
ax = axs[1] 

ngammas = 99 
av_r = Vector{Float64}(undef, ngammas)
av_sus = Vector{Float64}(undef, ngammas)
var_r = Vector{Float64}(undef, ngammas)
gammas = Vector{Float64}(undef, ngammas)
nsims = 100 

#Get the colormap
colors = ArtsyPalettes.met_brew("Isfahan1")
color = colors[5]

#Plot the mesoscopic data
data_meso     = "../../../data/diagrams/30harms_amplitude"


read_data_simpler!(data_meso, ngammas, nsims, av_r, av_sus, var_r, gammas)
lines!(ax, gammas, av_sus, color=color, label="Eqs. (8)")

#Plot the Kuramoto simulations...
ngammas = 43 
av_r = Vector{Float64}(undef, ngammas)
av_sus = Vector{Float64}(undef, ngammas)
var_r = Vector{Float64}(undef, ngammas)
gammas = Vector{Float64}(undef, ngammas)
nsims = 100 
data_kuramoto = "../../../data/diagrams/kuramoto_julia_bien"

read_data_simpler!(data_kuramoto, ngammas, nsims, av_r, av_sus, var_r, gammas)
scatter!(ax, gammas, av_sus, color=:black, label="Simulation")

ax.yticks = [0.0, 4e-4] 

xlims!(ax, 0.075, 0.125)
axislegend(ax, position=(0.9, 0.9))

ax.xlabel = "J"
ax.ylabel = "χ/N"



ax = axs[2] 

#Simulation data
ngammas = 100
av_r = Vector{Float64}(undef, ngammas)
av_sus = Vector{Float64}(undef, ngammas)
var_r = Vector{Float64}(undef, ngammas)
gammas = Vector{Float64}(undef, ngammas)

nsims = 3000 
data_path = "../../../data/gamma/offcritical_n1e6_sigsq0.1"
s2 = 0.1


read_data_simpler!(data_path, ngammas, nsims, av_r, av_sus, var_r, gammas)

#Get the colormap
colors = ArtsyPalettes.met_brew("Egypt")
colors = [colors[i] for i in [2,3]]

plot_effective_exponent(ax, gammas, av_sus)
axislegend(ax, position=(0.05, 1.5), padding=(0,0, 0,0))

#Add axes labels
xlabel = 0.94
ylabel = 0.75
StyleFuncs.label_axes(axs, pos=[xlabel ylabel; xlabel ylabel])

#Save the figure
save("figure2.pdf", fig, pt_per_unit = 1)