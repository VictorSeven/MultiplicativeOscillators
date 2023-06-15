module KuramotoCartesian
export get_timeseries, phase_diagram

using LinearAlgebra, Distributions

function get_corrs!(nharm, a, x, y, corr)
    for i=1:nharm
        for j=1:i-1
            if i+j <= nharm 
                xterm = x[i+j] 
                yterm = y[i+j]
            else
                xterm = 0.0 
                yterm = 0.0
            end

            corr[j,i] = j*i*a*(x[i-j] - xterm)
            corr[j+nharm,i] = j*i*a*(y[i-j] - yterm)
            corr[j+nharm,i+nharm] =  j*i*a*(x[i-j] + xterm)
        end

        if 2*i <= nharm
            corr[i,i] = i*i*a*(1.0 - x[2*i])
            corr[i+nharm,i] =  -i*i*a*y[2*i]
            corr[i+nharm,i+nharm] =  i*i*a*(1.0 + x[2*i])
        else
            corr[i,i] = i*i*a
            corr[i+nharm,i] = 0.0 
            corr[i+nharm,i+nharm] =  i*i*a
        end

    end
    corr .= Symmetric(corr, :L)
end

function get_correlated_vars!(corr, xi; thres=1e-10, maxits=10)

    good_corrs = isposdef!(corr)
    nits = 0
    while nits < maxits && !good_corrs 
        eigsys = eigen(corr)
        corr .= eigsys.vectors * (max.(eigsys.values, thres) .* eigsys.vectors')
        corr .= 0.5 * (corr + corr')
        eigsys = eigen(corr)
        good_corrs = !any(eigsys.values .< thres)
        nits += 1
    end

    #Now corr contains the Cholesky matrix,so 
    randn!(xi)
    xi .= corr * xi
    return nothing
end


function step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
    
    get_corrs!(nharm, a, old_x, old_y, corr)
    get_correlated_vars!(corr, xi)

    #Update k=1 so we can manually set old_r[0]=1 always
    detx = -0.5*(s2*old_x[1] + 2*w*old_y[1] + q*old_x[1]*(old_x[2] - 1.0) + q*old_y[1]*old_y[2])
    dety = 0.5*(-s2*old_y[1] + 2*w*old_x[1] + q*old_y[1]*(old_x[2] + 1.0) - q*old_x[1]*old_y[2])

    x[1] = old_x[1] + dt * detx + sqdt * xi[1] 
    y[1] = old_y[1] + dt * dety + sqdt * xi[1+nharm]

    #Update all k from 2 to nharm-1
    @simd for k=2:nharm-1
        detx = -0.5*k*(k*s2*old_x[k] + 2*w*old_y[k] + q*old_x[1]*(old_x[k+1] - old_x[k-1]) + q*old_y[1]*(old_y[k-1] + old_y[k+1]))
        dety = 0.5*k*(-k*s2*old_y[k] + 2*w*old_x[k] + q*old_y[1]*(old_x[k+1] + old_x[k-1]) + q*old_x[1]*(old_y[k-1] - old_y[k+1]))

        x[k] = old_x[k] + dt * detx + sqdt * xi[k]       
        y[k] = old_y[k] + dt * dety + sqdt * xi[k+nharm]
    end

    #Update last one to manually set old_r[nharm+1] = 0 
    k = nharm
    detx = -0.5*k*(k*s2*old_x[k] + 2*w*old_y[k] + q*old_x[1]*(0.0 - old_x[k-1]) + q*old_y[1]*(old_y[k-1]))
    dety = 0.5*k*(-k*s2*old_y[k] + 2*w*old_x[k] + q*old_y[1]*(0.0 + old_x[k-1]) + q*old_x[1]*(old_y[k-1]))

    x[k] = old_x[k] + dt * detx + sqdt * xi[k]       
    y[k] = old_y[k] + dt * dety + sqdt * xi[k+nharm]
end

function initial_conditions!(nharm, x, y)
    r = Vector{Float64}(undef, nharm)
    psi = Vector{Float64}(undef, nharm)

    r[1] = 0.01*rand()
    psi[1] = 2*π*rand() 
    k = collect(2:nharm)
    r[2:nharm] = r[1] .^ k 
    psi[2:nharm] = psi[1] .* k

    x .= r .* cos.(psi) 
    y .= r .* sin.(psi) 
end

function get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, fpath; w=0.01, dt=0.01, nsample=10)
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 

    old_x = Vector{Float64}(undef, nharm)
    old_y = Vector{Float64}(undef, nharm)
    x = Vector{Float64}(undef, nharm)
    y = Vector{Float64}(undef, nharm)

    corr = zeros(2*nharm, 2*nharm)

    initial_conditions!(nharm, old_x, old_y)

    for t=0:dt:t_thermal 
        step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
        old_x, old_y, x, y = x, y, old_x, old_y
    end

    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
            old_x, old_y, x, y = x, y, old_x, old_y

            if (nt % nsample == 0)
                write(output, "$t ")
                for k=1:nharm-1
                    rk = sqrt(x[k]^2 + y[k]^2)
                    write(output, "$rk ")
                end
                rk = sqrt(x[nharm]^2 + y[nharm]^2)
                write(output, "$rk\n")
            end
            t += dt
            nt += 1
        end
    end
    return abs(rk[1]) 
end

function phase_diagram(nharm, t_thermal, tf, q0, qf, nq, sys_size, s2, fpath; sampling=100, w=0.1, dt=0.01)
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 

    q_values = LinRange(q0, qf, nq)


    open(fpath, "w") do output 
        old_x = Vector{Float64}(undef, nharm)
        old_y = Vector{Float64}(undef, nharm)
        x = Vector{Float64}(undef, nharm)
        y = Vector{Float64}(undef, nharm)

        xi = zeros(2*nharm)
        for q in q_values


            corr = zeros(2*nharm, 2*nharm)

            initial_conditions!(nharm, old_x, old_y)

            for t=0:dt:t_thermal 
                step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
                old_x, old_y, x, y = x, y, old_x, old_y
            end

            avr = 0.0
            avr2 = 0.0
            nmeasures = 0

            t = 0
            nt = 0
            while t < tf 
                step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)

                if (nt % sampling == 0)
                    r = sqrt(x[1]^2 + y[1]^2) 
                    avr += r 
                    avr2 += r*r
                    nmeasures += 1
                end

                old_x, old_y, x, y = x, y, old_x, old_y


                t += dt
                nt += 1
            end

            avr /= nmeasures
            avr2 /= nmeasures

            write(output, "$q $avr $(avr2-avr*avr)\n")
        end
    end

    return nothing 
end


end #Module end