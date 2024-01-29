using CairoMakie
using DelimitedFiles

include("style_funcs.jl")
include("../core/theory_formulas.jl")

using .TheoryFormulas

function plot_thermodynamic_limit(axis)
    nparts = 4
    ngammas = 40 
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
            nextdata = readdlm("$data_path/diagrams_$(n)_verylong/diagram_$(index-1)_$(part-1)")
            #nextdata = readdlm("$data_path/diagrams_$(n)/diagram_$(index-1)_$(part-1)")
            data = vcat(data, nextdata) 
        end
        av_r += data[1:ngammas, 2] 
        av_sus += data[1:ngammas, 3] 
        gammas = data[1:ngammas, 1]
        
    end
    av_r /= nsims
    av_sus /= nsims
        
    teogammas = LinRange(0.0, 0.2, 100)
    r6  = r_6th(teogammas, s2)
    r2c = r_2th_cumulant(teogammas, s2)
    r2c[teogammas .< 0.1] .= 0

    r   = r_oa(teogammas, s2)
    datatyul = readdlm("dtyul")
    rtrue = readdlm("../../../data/gamma/theoretical/sf_50harms")[:,1] #integrate_full.(0.1, gammas, s2)

    lines!(axis, teogammas, r6, label="6th harmonic")
    lines!(axis, teogammas, r2c, label="3rd cumulant")
    lines!(axis, datatyul[:,1], datatyul[:,2], label="Tyulkina et al.")
    lines!(axis, teogammas, rtrue, label="Full System")
    lines!(axis, teogammas, r, label="Ott-Antonsen")

    scatter!(axis, gammas, av_r, markersize=4, color=:gray, label="Simulation")

    #axislegend(axis, position=(0.01, 0.5), framevisible=false, rowgap=0, patchlabelgap=1, patchsize=(10,10)) #linepoints=[Point2f(0.5, 0.5), Point2f(1.0, 0.5)] )
    create_legend(axis, (0.005, 0.9))

    axis.xlabel = L"J"
    axis.ylabel = L"R"

end

function plot_finite_sizes(axis)
    nparts = 3
    ngammas = 50
    nsims = 100
    data_path = "../../../data/diagrams"

    s2 = 0.1

    av_r = Vector{Float64}(undef, ngammas)
    av_sus = Vector{Float64}(undef, ngammas)
    gammas = Vector{Float64}(undef, ngammas)

    system_sizes = [100, 1000, 10000]
    cmap =  cgrad(:davos)
    colors = [cmap.colors[i] for i in [20,50,100]]

    i=1
    for n in system_sizes
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
        scatter!(axis, gammas, av_r, markersize=4, color=colors[i], label=L"N=%$n")
        r_th = TheoryFormulas.finite_size_r(gammas, s2, n)
        lines!(axis, gammas, r_th, color=colors[i])
        i += 1
    end

    axis.xlabel = L"J"

    create_legend(axis, (0.005, 0.9))
end


#lines!(data[:,1], av_r)

fig = Figure(resolution=two_col_size(2*1.618), fontsize=9, figure_padding=7)
group = fig[1,1] = GridLayout()
axs = [Axis(group[1,j], xgridvisible=false, ygridvisible=false) for j=1:2]
#, xticks=[0.0, 0.1, 0.2]

colgap!(group, 10)
hidespines!.(axs, :t, :r)

plot_thermodynamic_limit(axs[1])
plot_finite_sizes(axs[2])

hideydecorations!(axs[2])
for ax in axs
    ax.xticks = [0.0, 0.1, 0.2] 
    ylims!(ax, 0.0, 0.85)
end
#text!(axs[1], 0.1,0.5; text="hola")
label_axes(axs)
save("figure1.pdf", fig, pt_per_unit = 1)