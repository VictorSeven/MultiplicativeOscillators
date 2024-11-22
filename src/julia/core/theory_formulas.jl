module TheoryFormulas

#We will need this
using Statistics

#Export some functions to be used
export r_6th, r_2th_cumulant, r_oa, finite_size_r
export integrate_tyulkina, diagram_tyulkina
export integrate_full, integrate_amplitudes, diagram_amplitudes


# ----- Simple, direct math formulas ---- 

"""
    r_oa(q::Float64, s2::Float64)

The theoretical form of the Kuramoto order parameter under
the Ott-Antonsen ansatz given a coupling q and noise intensity squared s2
"""
function r_oa(q, s2)
    return @. real(sqrt(complex((q-s2)/q))) 
end

"""
    r_6th(q::Float64, s2::Float64)

Theoretical formula for the Kuramoto model up to 6th harmonic
given a coupling q and noise intensity squared s2
"""
function r_6th(q, s2)
    ratio = @. s2 * sqrt(2) / q
    root = @. sqrt(complex(51*q*q - 132 * q * s2 + 306*s2*s2)) 
    numerator = @. complex(24*s2 - 9*q - root)
    denominator = @. complex(q - 9*s2)
    return @. real(ratio * sqrt(numerator / denominator))
end

"""
    r_2th_cumulant(q::Float64, s2::Float64)

Theoretical formula for the Kuramoto model by cancelling the third 
circular cumulant 
"""
function r_2th_cumulant(q, s2)
    root = @. sqrt(complex(4*q*q + 4*q*s2 - 7*s2*s2))
    numerator = @. complex(2*q - 3*s2 + root) 
    return @. real(0.5 * sqrt(numerator / q))
end

"""
    finite_size_r(q::Float64, s2::Float64, n::Int)

Theoretical solutiion for the Ott-Antonsen ansatz 
using a finite number of oscillators, N
"""
function finite_size_r(q, s2, n)
    aux = @. q * n - s2 * (n - 1)  
    aux = @. aux + sqrt(4*n*q*s2 + aux^2)
    return @. sqrt(aux / (2*q*n))
end


# ---- Obtained by the integration of an ODE system ---- 

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
    amplitude_finitesize(nharm::Integer, sys_size::Integer, q::Float64, s2::Float64, oldr::Float64, r::Float64, dt::Float64)

Definition of the dynamical system up to `nharm` harmonics, assuming a system 
size `sys_size`, a coupling `q` and noise intensity square `s2`.
The method also takes the old values `oldz` and rewrites the 
new `z` ones. It takes the timestep `dt`.
"""
function amplitude_finitesize!(nharm::Integer, sys_size::Integer, q::Float64, s2::Float64, oldr::AbstractVector, r::AbstractVector, dt::Float64)
    k = 1 
    det = 0.5*k*(q*oldr[1]*(1.0 - oldr[k+1]) -k*s2*oldr[k]) + k^2*s2*(1 + oldr[2k]) / (4 * sys_size * oldr[1])
    r[k] = oldr[k] + dt * det 

    #From harmonic 2 to nharm-1
    @simd for k=2:nharm-1
        det = 0.5*k*(q*oldr[1]*(oldr[k-1] - oldr[k+1]) -k*s2*oldr[k]) 
        r[k] = oldr[k] + dt * det 
    end

    #Update the last harmonic by hand to include the closure condition z[nharm+1]=0
    k = nharm 
    det = 0.5*k*(q*oldr[1]*oldr[k-1] - k*s2*oldr[k]) 
    r[k] = oldr[k] + dt * det 
end

"""
    integrate_amplitudes(nharms::Integer, sys_size::Integer, q::Float64, s2::Float64; dt=0.01, tf=500.0)

Simple integration of the theory for `nharms` harmonics, 
a coupling `q` and noise intensity square `s2`. Optionally, one can
also set the timestep `dt`, and the simulation time `tf`
"""
function integrate_amplitudes(nharms::Integer, sys_size::Integer, q::Float64, s2::Float64; dt=0.01, tf=500.0)
    #Initialize the angles and get Kuramoto order parameter 
    #and its fluctuations
    oldr = 0.1*rand(nharms)
    r = zeros(nharms)

    #Then just integrate everything for the specified time
    nits = Int(tf/dt)
    for i=1:nits
        amplitude_finitesize!(nharms, sys_size, q, s2, oldr, r, dt) 
        oldr, r = r, oldr
    end

    #Return the result
    return oldr 
end

"""
    diagram_amplitudes(nharms::Integer, sys_size::Integer, q0::Float64, qf::Float64, nq::Integer, s2::Float64, path::String; tf=5000.0)

Generates a phase diagram using the new theory including FS corrections at size `sys_size`.
From coupling `q0` to `qf` making `nq` subdivisions, and square noise intensity `s2`. 
Results are then stored in `path`. One can also specify the time for the integration, 
We set `tf=5000.0` to ensure thermalization
"""
function diagram_amplitudes(nharms::Integer, sys_size::Integer, q0::Float64, qf::Float64, nq::Integer, s2::Float64, path::String; tf=5000.0, dt=0.01, writetofile=true)
    #Open the file and iterate over the couplings
    if writetofile 
        open(path, "w") do output
            for q in LinRange(q0, qf, nq)
                r = integrate_amplitudes(nharms, sys_size, q, s2; dt=dt, tf=tf)
                write(output, "$q $(r[1])\n")
            end
        end
    else
        rvalues = zeros(nq)
        i = 1
        qlist = LinRange(q0, qf, nq)
        for q in qlist 
            rk = integrate_amplitudes(nharms, sys_size, q, s2; dt=dt, tf=tf)
            rvalues[i] = rk[1]
            i += 1
        end
        return qlist, rvalues
    end
end


"""
    tyulkina(w::Float64, q::Float64, oldz::Float64, oldx::Float64, dt::Float64)

Definition of the dynamical system representing the Kuramoto order
parameter in Tyulkina et al. PRL 2018 (with non-zero frequency `w`),
for a coupling `q` and noise intensity square `s2`. 
The method takes also values at iteration `t-dt` (`oldz` and `oldx`)
as well as the timestep `dt`.
"""
function tyulkina(w, q, s2, oldz, oldx, dt)
    z = oldz + dt * ( 0.5*(q-s2+2*im*w)*oldz - 0.5*q*(oldz * abs2(oldz) + conj(oldz)*oldx)) 
    x = oldx + dt * ( 2*im*w*oldx - s2*(2*oldx + oldz^2) - 2*q*oldx*abs2(oldz) )
    return z,x
end

"""
    diagram_tyulkina(w::Float64, q::Float64, s2::Float64g; tf=500.0)    

Simple integration of the Tyulkina et al. ODEs for a frequency `w`, 
a coupling `q` and noise intensity square `s2`. Optionally, one can
also set the timestep `dt`, the number of oscillators `n` (to specify
initial density), and the simulation time `tf`
"""
function integrate_tyulkina(w, q, s2; dt=0.01, n=100000, tf=500.0)
    #Initialize the angles and get Kuramoto order parameter 
    #and its fluctuations
    angles = 2π*rand(n)
    oldz = mean(exp.(im*angles))
    z2 = mean(exp.(2*im*angles))
    oldx = z2 - oldz.^2  

    #Then just integrate everything for the specified time
    nits = Int(tf/dt)
    z = 0.0+0.0im 
    for i=1:nits
        z, x = tyulkina(w, q, s2, oldz, oldx, dt)
        oldz, oldx = z, x
    end

    #Return the result
    return abs(oldz)
end


"""
    diagram_tyulkina(q0::Float64, qf::Float64, nq::Int, w::Float64, s2::Float64, path::String; tf=5000.0)    

Generates a phase diagram using the system by Tyulkina et al.
from coupling `q0` to `qf` making `nq` subdivisions, given a frequency
`w` and square noise intensity `s2`. Results are then stored in `path`.
One can also specify the time for the integration, 
we set `tf=5000.0` to ensure thermalization
"""
function diagram_tyulkina(qspace, w, s2; tf=5000.0)
    r = Vector{Float64}(undef, length(qspace))
    for (i,q) in enumerate(qspace) 
        r[i] = integrate_tyulkina(w, q, s2; tf=tf)
    end
    return r
end



"""
    timeseries_tyulkina(w::Float64, q::Float64, s2::Float64, path::String; dt=0.01, tf=1000.0)    

Generates a timeseries of the system by Tyulkina et al., given a frequency
`w` and square noise intensity `s2`. Results are then stored in `path`.
One can also specify the time for the integration, 
"""
function timeseries_tyulkina(w, q, s2, outputpath; dt=0.01, tf=1000.0)
    #Initialize random angles
    angles = 2π*rand(100000)
    oldz = mean(exp.(im*angles))
    z2 = mean(exp.(2*im*angles))
    oldx = z2 - oldz.^2  

    #Do the integration
    z = 0.0
    nits = 0
    open(outputpath, "w") do output
        for t=0.0:dt:tf 
            z, x = tyulkina(w, q, s2, oldz, oldx, dt)
            oldz, oldx = z, x
            #Do not write all the results, just 1 each 10
            if nits%10==0
                write(output, "$t $(abs(oldz))\n")
            end
            nits += 1
        end
    end
end


end