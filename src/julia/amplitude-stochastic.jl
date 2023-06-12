using LinearAlgebra
using Distributions

function get_corrs!(nharm, a, r, corr)
    for i=1:nharm
        for j=1:i-1
            if i+j < nharm 
                term = r[i+j]
            else
                term = 0.0 
            end

            corr[j,i] = i*j*a*(r[i-j] - term)  
        end

        if 2*i < nharm
            corr[i,i] = i*i*a*(1.0 - r[2*i])
        else
            corr[i,i] = i*i*a
        end

    end
    corr .= Symmetric(corr, :L)
end

function ensure_corr!(corr, thres=1e-10, maxits=10)
    good_corrs = isposdef(corr)
    nits = 0
    while nits < maxits && !good_corrs 
        eigsys = eigen(corr)
        corr .= eigsys.vectors * (max.(eigsys.values, thres) .* eigsys.vectors')
        corr .= 0.5 * (corr + corr')
        eigsys = eigen(corr)
        good_corrs = !any(eigsys.values .< thres)
        nits += 1
    end
end


function step!(nharm, oldr, r, corr, w, q, s2, a, dt, sqdt, t)
    
    get_corrs!(nharm, a, oldr, corr)
    ensure_corr!(corr)
    mvn = MvNormal(zeros(nharm), corr) 
    xi = rand(mvn, 1)

    k = 1 #z[0]=1
    det = 0.5*k*(q*oldr[1]*(1.0 - oldr[k+1]) -k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    #r[k] = max(r[k], 0.0)
    @simd for k=2:nharm-1
        det = 0.5*k*(q*oldr[1]*(oldr[k-1] - oldr[k+1]) -k*s2*oldr[k]) 
        r[k] = oldr[k] + dt * det + sqdt * xi[k]
        #r[k] = max(r[k], 0.0)
    end
    k = nharm #z[nharm+1]=0
    det = 0.5*k*(q*oldr[1]*oldr[k-1] - k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    #r[k] = max(r[k], 0.0)
end

function initial_conditions!(nharm)
    r = Vector{Float64}(undef, nharm)

    r[1] = 0.01*rand()
    k = collect(2:nharm)
    r[2:nharm] = r[1] .^ k 

    return r
end

function get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, fpath, detfs=true)
    a = 0.5 * s2 / sys_size
    dt = 0.01 
    w = 0.01
    sqdt = sqrt(dt) 
    sampling = 0.1 
    nsample = floor(sampling/dt)

    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)

    corr = zeros(2*nharm, 2*nharm)

    old_r = initial_conditions!(nharm)

    for t=0:dt:t_thermal 
        step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t)
        old_r, r = r, old_r
    end

    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t)
            old_r, r = r, old_r

            if (nt % nsample == 0)
                write(output, "$t ")
                for k=1:nharm-1
                    write(output, "$(r[k]) ")
                end
                write(output, "$(r[k])\n")
            end
            t += dt
            nt += 1
        end
    end
    return abs(rk[1]) 
end

function compute_stat_r(nharm, trelax, tf, q, sys_size, s2)
    a = 0.5 * s2 / sys_size
    dt = 0.01 
    w = 0.01
    sqdt = sqrt(dt) 

    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)


    corr = zeros(nharm, nharm)

    old_r = initial_conditions!(nharm)

    for t=0:dt:trelax
        step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t)
        old_r, r = r, old_r
    end

    avr = 0.0
    avr2 = 0.0
    for t=0:dt:tf
        step!(nharm, old_r, r, corr, w, q, s2, a, dt, sqdt, t)
        avr += r[1]
        avr2 += r[1]*r[1]

        old_r, r = r, old_r
    end

    avr /= tf/dt
    avr2 /= tf/dt

    return abs(avr), avr2 - avr*avr 
end