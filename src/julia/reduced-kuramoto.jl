using LinearAlgebra
using Distributions

function get_corrs!(nharm, a, r, corr)
    for i=1:nharm
        for j=1:i-1
            if i+j < nharm
                ipj = r[i+j]
            else
                ipj = 0.0
            end
            corr[j,i] = i*j*a*(r[i-j] - ipj)
        end 

        if 2*i<nharm
            corr[i,i] = i*i*a*(1.0 - r[2*i])
        else 
            corr[i,i] = i*i*a
        end
    end
    corr .= Symmetric(corr)
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


function steptl!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
    
    get_corrs!(nharm, a, old_r, corr)
    ensure_corr!(corr)
    mvn = MvNormal(zeros(nharm), corr) 
    xi = rand(mvn, 1)


    #Update k=1 so we can manually set old_r[0]=1 always
    det = 0.5*(q * old_r[1] * (1.0 - old_r[2]) - s2*old_r[1]) 
    r[1] = old_r[1] + dt * det + sqdt * xi[1]
    r[1] = max(r[1], 0)

    #Update all k from 2 to nharm-1
    for k=2:nharm-1
        det = 0.5*k* (q * old_r[1] * (old_r[k-1] - old_r[k+1]) - k*s2*old_r[k]) 
        r[k] = old_r[k] + dt * det + sqdt * xi[k]
        r[k] = max(r[k], 0)
    end

    #Update last one to manually set old_r[nharm+1] = 0 
    k = nharm
    det = 0.5*k* (q * old_r[1] * old_r[k-1] - k*s2*old_r[k]) 
    r[k] = old_r[k] + dt * det + sqdt * xi[k]
    r[k] = max(r[k], 0)
end

function stepfs!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
    
    get_corrs!(nharm, a, old_r, corr)
    ensure_corr!(corr)
    mvn = MvNormal(zeros(nharm), corr) 
    xi = rand(mvn, 1)


    #Update k=1 so we can manually set old_r[0]=1 always
    det = 0.5*(q * old_r[1] * (1.0 - old_r[2]) - s2*old_r[1]) + a * (1+old_r[1]^2)/old_r[1] 
    r[1] = old_r[1] + dt * det + sqdt * xi[1]
    r[1] = max(r[1], 1e-20)
    println("1 ", old_r[1], " ", r[1], " ", det)

    #Update all k from 2 to nharm-1
    for k=2:nharm-1
        det = 0.5*k* (q * old_r[1] * (old_r[k-1] - old_r[k+1]) - k*s2*old_r[k]) + a * (1+old_r[k]^2)/old_r[k] 
        r[k] = old_r[k] + dt * det + sqdt * xi[k]
        r[k] = max(r[k], 1e-20)
        println(k, " ", old_r[k], " ", r[k], " ", det)
    end

    #Update last one to manually set old_r[nharm+1] = 0 
    k = nharm
    det = 0.5*k* (q * old_r[1] * old_r[k-1] - k*s2*old_r[k]) + a * (1+old_r[k]^2)/old_r[k] 
    r[k] = old_r[k] + dt * det + sqdt * xi[k]
    r[k] = max(r[k], 1e-20)
    println(k, " ", old_r[k], " ", r[k], " ", det)
    println()
end

function initial_conditions!(nharm, r)
    r[1] = 0.2 + 0.5 * rand()
    r[2:nharm] = r[1] .^ collect(2:nharm)
end

function get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, fpath, detfs=true)
    a = 0.5 * s2 / sys_size
    dt = 0.01 
    sqdt = sqrt(dt) 
    sampling = 1e-2
    nsample = floor(dt/sampling)


    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)

    corr = zeros(nharm, nharm)

    initial_conditions!(nharm, old_r)

    stepfun! = detfs ? stepfs! : steptl!

    for t=0:dt:t_thermal 
        #steptl!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
        stepfun!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
        old_r, r = r, old_r
    end

    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            #steptl!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
            stepfun!(nharm, old_r, r, corr, q, s2, a, dt, sqdt, t)
            old_r, r = r, old_r

            if (nt % nsample == 0)
                write(output, "$t ")
                for k=1:nharm-1
                    rk = r[k]
                    write(output, "$rk ")
                end
                rk = r[nharm]
                write(output, "$rk\n")
            end
            t += dt
            nt += 1
        end
    end
    return corr 
end