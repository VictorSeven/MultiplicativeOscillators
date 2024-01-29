module AmplitudeEquations
export get_timeseries, phase_diagram 

using LinearAlgebra, Distributions, Random

function get_corrs!(nharm, a, r, corr)
    for i=1:nharm
        for j=1:i-1
            if i+j <= nharm 
                term = r[i+j]
            else
                term = 0.0 
            end

            corr[j,i] = i*j*a*(r[i-j] - term)  
        end

        if 2*i <= nharm
            corr[i,i] = i*i*a*(1.0 - r[2*i])
        else
            corr[i,i] = i*i*a
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


function step!(nharm, oldr, r, corr, w, q, s2, a, dt, sqdt, t, xi)
    
    get_corrs!(nharm, a, oldr, corr)
    get_correlated_vars!(corr, xi)

    k = 1 #z[0]=1
    det = 0.5*k*(q*oldr[1]*(1.0 - oldr[k+1]) -k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    r[k] = max(r[k], 0.0)
    @simd for k=2:nharm-1
        det = 0.5*k*(q*oldr[1]*(oldr[k-1] - oldr[k+1]) -k*s2*oldr[k]) 
        r[k] = oldr[k] + dt * det + sqdt * xi[k]
        r[k] = max(r[k], 0.0)
    end
    k = nharm #z[nharm+1]=0
    det = 0.5*k*(q*oldr[1]*oldr[k-1] - k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    r[k] = max(r[k], 0.0)
end

function initial_conditions!(nharm, r)
    r[1] = 0.01*rand()
    k = collect(2:nharm)
    r[2:nharm] = r[1] .^ k 

    return nothing
end

function get_timeseries(nharm, t_thermal, tf, w, q, sys_size, s2, fpath; dt=0.01, nsample=10)
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 
    xi = zeros(nharm)

    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)

    corr = zeros(nharm, nharm)

    initial_conditions!(nharm, old_r)

    for t=0:dt:t_thermal 
        step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t,  xi)
        old_r, r = r, old_r
    end

    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t, xi)
            old_r, r = r, old_r

            if (nt % nsample == 0)
                write(output, "$t ")
                for k=1:nharm-1
                    write(output, "$(r[k]) ")
                end
                write(output, "$(r[nharm])\n")
            end
            t += dt
            nt += 1
        end
    end
    return abs(r[1]) 
end


function phase_diagram(nharm, t_thermal, tf, q0, qf, nq, sys_size, s2, fpath; sampling=100, w=0.01, dt=0.01)
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 

    open(fpath, "w") do output 
        
        q_values = LinRange(q0, qf, nq)
        xi = zeros(nharm)

        old_r = Vector{Float64}(undef, nharm)
        r = Vector{Float64}(undef, nharm)

        for q in q_values

            corr = zeros(nharm, nharm)

            initial_conditions!(nharm, old_r)

            for t=0:dt:t_thermal 
                step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t, xi)
                old_r, r = r, old_r
            end

            avr = 0.0
            avr2 = 0.0
            nmeasures = 0

            t = 0
            nt = 0
            while t < tf 
                step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t, xi)

                if (nt % sampling == 0)
                    avr += r[1] 
                    avr2 += r[1]*r[1]
                    nmeasures += 1
                end

                old_r, r = r, old_r


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