"""
Contains code to integrate the amplitude equations of the paper, 
i.e., after assuming the ansatz for the angles.
"""
module FixedPointNoise 
export get_timeseries, phase_diagram 

using LinearAlgebra, Distributions, Random, NearestCorrelationMatrix


"""
    step!(nharm::Int, oldr::Vector, r::Vector, corr:Matrix, w::Float64, q::Float64, s2::Float64, a::Float64, dt::Float64, sqdt::Float64, t::Float64, xi::Vector)

An integration step using `nharm` harmonics, computing the new Kuramoto-Daido parameters `r` from
the old ones `oldr`. Uses the frequency `w`, coupling `q`, noise intensity `s2` 
square and mesoscopic noise `a`. The variables `corr` and `xi` are the correlation matrix and stochastic
noise and are given to avoid allocating memory. The method also needs the timestep `dt`, its 
square root `sqdt` and current timestep `t`.
"""
function step!(nharm, oldr, r, w, q, s2, g, dt, sqdt, t, xi)
    #Directly generate diagonal multiplicative noise 
    randn!(xi)
    xi .= g * xi

    #Computation of the first harmonic 
    #We include 0th harmonic manually, since the 0th is always r[0]=1)
    #Then clamp the order parameter to its validity range
    k = 1 
    det = 0.5*k*(q*oldr[1]*(1.0 - oldr[k+1]) -k*s2*oldr[k])  
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    r[k] = min(1.0, max(r[k], 0.0))

    #From harmonic 2 to nharm-1
    @simd for k=2:nharm-1
        det = 0.5*k*(q*oldr[1]*(oldr[k-1] - oldr[k+1]) -k*s2*oldr[k]) 
        r[k] = oldr[k] + dt * det + sqdt * xi[k] 
        r[k] = min(1.0, max(r[k], 0.0))
    end

    #Update the last harmonic by hand to include the closure condition z[nharm+1]=0
    k = nharm 
    det = 0.5*k*(q*oldr[1]*oldr[k-1] - k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det + sqdt * xi[k]
    r[k] = min(1.0, max(r[k], 0.0))
end

function construct_gmatrix(nharm, system_size, w, q, s2, dtsys=0.001)

    r0 = rand(2*nharm)
    psi0 = 2π * rand(2*nharm)
    oldz = r0 .* exp.(1im*psi0) 
    #z = Vector{ComplexF64}(undef, 2*nharm) 
    z = zeros(2*nharm) + im*zeros(2*nharm)

    tf = 500.0
    nsteps = floor(tf / dtsys)

    for i=1:nsteps
        full_system(2*nharm, w, q, s2, oldz, z, dtsys)
        oldz, z = z, oldz
    end

    r = abs.(oldz)
   
    corr = Matrix{Float64}(undef, nharm, nharm)
    d = zeros(nharm)
    get_corrs!(nharm, r, corr, d)

    cholesky!(corr)
    tril!(corr)
    
    corr .=  sqrt(0.5*s2/system_size) * Diagonal(d) * corr 
    
    return r,  corr  
end

"""
    full_system(nharm::Int, w::Float64, q::Float64, s2::Float64, oldz::Float64, z::Float64, dt::Float64)

Definition of the dynamical system up to `nharm` harmonics, for 
frequency `w`, a coupling `q` and noise intensity square `s2`.
The method also takes the old values `oldz` and rewrites the 
new `z` ones. It takes the timestep `dt`.
"""
function full_system(nharm, w, q, s2, oldz, z, dt)
    #First harmonic is k=1 since z[0]=1 by definition
    #Do it outside of the loop to add this initial condition
    k = 1 
    z[k] = oldz[k] + dt *( oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*(oldz[1] - conj(oldz[1])*oldz[k+1]))
    
    #Loop over everyone except the last...
    @simd for k=2:nharm-1
        z[k] = oldz[k] + dt *( oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*(oldz[1]*oldz[k-1] - conj(oldz[1])*oldz[k+1]))
    end

    #So we can manually add the closure condition z[nharm+1]=0 for the last
    k = nharm 
    z[k] = oldz[k ] + dt * (oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*oldz[1]*oldz[k-1] )
end

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
            #Matrix element
            corr[j,i] = i*j*(r[i-j] - r[i+j])  
            corr[i,j] = corr[j,i] 
        end

        #Fill the diagonal
        corr[i,i] = i*i*(1.0 - r[2*i]) 
        d[i] = sqrt(corr[i,i])
    end


    invd = 1. ./ d

    #Enforce symmetry
    corr .= Diagonal(invd)*corr*Diagonal(invd) 
    corr .= Symmetric(corr)
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
    sqdt = sqrt(dt) 

    #Initialize vectors
    r = Vector{Float64}(undef, nharm)
    xi = zeros(nharm)

    #Evaluate the correlation matrix at the stationary solution and decompose it to get the g
    old_r, g = construct_gmatrix(nharm, sys_size, w, q, s2)

    #Write results to the file up to time tf
    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_r, r, w, q, s2, g, dt, sqdt, t, xi)
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


end #Module end