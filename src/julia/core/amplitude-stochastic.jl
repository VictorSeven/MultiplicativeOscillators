"""
Contains code to integrate the amplitude equations of the paper, 
i.e., after assuming the ansatz for the angles.
"""
module AmplitudeEquations
export get_timeseries, phase_diagram 

using LinearAlgebra, Distributions, Random, NearestCorrelationMatrix

"""
    get_corrs!(nharm::Int, a::Float64, r::Vector, corr::Matrix)

Computes the multiplicative noise matrix and stores it in `corr`, assuming
`nharm` harmonics, noise intensity `a` and Kuramoto-Daido parameters `r`
"""
function get_corrs!(nharm, r, corr, d)
    #d = Vector{Float64}(undef, nharm)
    
    #Fill half of the matrix (a corr matrix is symmetric)
    for i=1:nharm
        for j=1:i-1
            #Closure condition. All harmonics larger than nharm are 0
            term = i+j<=nharm ? r[i+j] : r[1]^(i+j) 

            #Matrix element
            corr[j,i] = i*j*(r[i-j] - term)  
        end

        #Fill the diagonal
        corr[i,i] = 2*i <= nharm ? i*i*(1.0 - r[2*i]) : i*i*(1.0 - r[1]^(2*i))
        #d[i] = 1/sqrt(corr[i,i])
        d[i] = sqrt(corr[i,i])
    end

    invd = 1. ./ d

    #Enforce symmetry
    corr .= Diagonal(invd)*corr*Diagonal(invd) 
    corr .= Symmetric(corr)
end

"""
    get_correlated_vars!(corr::Matrix, xi::Vector; thres=1e-10, maxits=10)

Obtain the stochastic noise from the correlation matrix `corr` and store it in the 
vector `xi`. This method will try to avoid ill-defined correlation matrices by capping
eigenvalues smaller than `thres` to `thres` (ensures positive definite) and symmetrizing 
after eigenvalues are capped. The check is repeated after symmetrization, for a total of
`maxits` iterations.
"""
function get_correlated_vars!(corr, sqa, t, xi, d)
    #Is our matrix positive definite? 
    good_corrs = isposdef(corr)

    #If it is, get our random variables
    if good_corrs
        #Cholesky decomposition. This fills corr with no zeros, so take just the lower triangle
        cholesky!(corr)
        tril!(corr)

        #Random variables
        randn!(xi)

        #Recover the original variance
        #xi .=  sqrt(a) * Diagonal(d) * corr * xi
        xi .=  sqa * d .* (corr * xi)
    else 
        #If not, get the nearest (normalized) correlation matrix and then do the same as above
        nearest_cor!(corr)
        cholesky!(corr)
        tril!(corr)
        randn!(xi)
        #xi .=  sqrt(a) * Diagonal(d) * corr * xi
        xi .=  sqa * d .* (corr * xi)
    end
    return nothing
end

"""
    step!(nharm::Int, oldr::Vector, r::Vector, corr:Matrix, w::Float64, q::Float64, s2::Float64, a::Float64, dt::Float64, sqdt::Float64, t::Float64, xi::Vector)

An integration step using `nharm` harmonics, computing the new Kuramoto-Daido parameters `r` from
the old ones `oldr`. Uses the frequency `w`, coupling `q`, noise intensity `s2` 
square and mesoscopic noise `a`. The variables `corr` and `xi` are the correlation matrix and stochastic
noise and are given to avoid allocating memory. The method also needs the timestep `dt`, its 
square root `sqdt` and current timestep `t`.
"""
function step!(nharm, oldr, r, corr, sys_size, q, s2, sqa, dt, sqdt, t, xi)
    #Get the correlation matrix and generate multiplicative noise
    d = Vector{Float64}(undef, nharm)
    #get_corrs!(nharm, a, oldr, corr)
    get_corrs!(nharm, oldr, corr, d)
    get_correlated_vars!(corr, sqa, t, xi, d)

    #Computation of the first harmonic 
    #We include 0th harmonic manually, since the 0th is always r[0]=1)
    #Then clamp the order parameter to its validity range
    #In order to take small systems into account, we include finite size correction
    #so it is not necessary to clamp r to be always higher than 0 (that causes problems)
    k = 1 
    det = 0.5*k*(q*oldr[1]*(1.0 - oldr[k+1]) -k*s2*oldr[k]) + k^2*s2*(1 + oldr[2k]) / (4 * sys_size * oldr[1])
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    r[k] = min(1.0, r[k])
    r[k] = abs(r[k])

    #From harmonic 2 to nharm-1. No finite-size correction from here. 
    @inbounds @simd for k=2:nharm-1
        det = 0.5*k*(q*oldr[1]*(oldr[k-1] - oldr[k+1]) -k*s2*oldr[k]) 
        r[k] = oldr[k] + dt * det + sqdt * xi[k]
        r[k] = min(1.0, max(r[k], 0.0))
    end

    #Update the last harmonic by hand to include the closure condition z[nharm+1]=0
    k = nharm 
    det = 0.5*k*(q*oldr[1]*oldr[k-1] - k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    #r[k] = max(r[k], 0.0)
    r[k] = min(1.0, max(r[k], 0.0))
end

"""
    initical_conditions!(nharm::Int, r::Vector)

Set the initial conditions, filling the vector `r` (of size `nharm`) of Kuramoto-Daido amplitudes.
"""
function initial_conditions!(nharm, r)
    r[1] = 0.01*rand()
    k = collect(2:nharm)
    r[2:nharm] = r[1] .^ k 

    return nothing
end

"""
    get_timeseries(nharm::Int, t_thermal::Float64, tf::Float64, w::Float64, q::Float64, sys_size::Int, s2::Float64, fpath::String; dt=0.01, nsample=10 )

Make a timeseries with a system of `nharm` harmonics, including frequency `w`, coupling `q`, system size `sys_size`
and noise intensity squared `s2`. The system will relax for a time `t_thermal` and then write the 
values of all the Kuramoto-Daido amplitudes to `fpath` for a time `tf`. One can set timestep `dt` and
the number of iterations between writing to the file `nsample`.

Output file format: time r[1] r[2] ... r[nharm]
"""
function get_timeseries(nharm, t_thermal, tf, w, q, sys_size, s2, fpath; dt=0.01, nsample=10)
    #Get constants we need for the integration
    sqa = sqrt(0.5 * s2 / sys_size)
    sqdt = sqrt(dt) 

    #Initialize vectors
    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)
    xi = zeros(nharm)
    corr = zeros(nharm, nharm)

    #Initial conditions
    initial_conditions!(nharm, old_r)

    #Thermalize
    for t=0:dt:t_thermal 
        step!(nharm, old_r, r, corr, sys_size, q, s2, sqa, dt, sqdt, t,  xi)
        old_r, r = r, old_r
    end

    #Write results to the file up to time tf
    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_r, r, corr, sys_size, q, s2, sqa, dt, sqdt, t, xi)
            old_r, r = r, old_r

            #Do not write all the iterations to avoid too dense file if dt is small
            if (nt % nsample == 0)
                write(output, "$t ")
                #Write all Kuramoto-Daido amplitudes by columns
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

"""
    phase_diagram(nharm::Int, t_thermal::Float64, tf::Float64, q0::Float64, qf::Float64, nq::Int, sys_size::Int, s2::Float64, fpath::String; sampling=100, w=0.01, dt=0.01)

Compute a phase diagram for the Kuramoto order parameter using  the system with `nharm` harmonics, 
leaving `t_thermal` time units to thermalize and measuring over `tf` time units. 
The phase diagram is done from a coupling `q0` to  a coupling `qf` in `nq` divisions. 
The system will have noise square intensity `s2` and system size `sys_size` and will write the 
results to `fpath`. The `sampling` variable gives the number of iterations between measurements 
to avoid correlations; system has a fixed frequency `w` and timestep `dt`.

Output file format: q mean_r var_r
"""
function phase_diagram(nharm, t_thermal, tf, q0, qf, nq, sys_size, s2, fpath; sampling=100, w=0.01, dt=0.01)
    #Constants of the system
    sqa = sqrt(0.5 * s2 / sys_size)
    sqdt = sqrt(dt) 

    #Initialize the vectors
    q_values = LinRange(q0, qf, nq)
    xi = zeros(nharm)

    old_r = Vector{Float64}(undef, nharm)
    r = Vector{Float64}(undef, nharm)
    
    corr = zeros(nharm, nharm)
    #Start simulations
    open(fpath, "w") do output 
        for q in q_values
            #Initial conditions
            initial_conditions!(nharm, old_r)

            #Thermalize
            for t=0:dt:t_thermal 
                step!(nharm, old_r, r, corr, sys_size, q, s2, sqa, dt, sqdt, t, xi)
                old_r, r = r, old_r
            end

            #Prepare to measure
            avr = 0.0
            avr2 = 0.0
            nmeasures = 0

            #Measure over tf time units
            t = 0
            nt = 0
            while t < tf 
                step!(nharm, old_r, r, corr, sys_size, q, s2, sqa, dt, sqdt, t, xi)

                #Here we take a measure, every sampling iterations
                if (nt % sampling == 0)
                    avr += r[1] 
                    avr2 += r[1]*r[1]
                    nmeasures += 1
                end

                #Prepare for next iteration
                old_r, r = r, old_r

                t += dt
                nt += 1
            end

            #Finish the average
            avr /= nmeasures
            avr2 /= nmeasures

            #Write result
            write(output, "$q $avr $(avr2-avr*avr)\n")
        end
    end

    return nothing 
end


end #Module end