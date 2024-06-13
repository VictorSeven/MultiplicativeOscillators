using Random 

"""
    initial_conditions!(N::Integer, phi::AbstractVector)

Sets random initial conditions for the phases and computes the initial Kuramoto order parameter.

    Return: R, psi
"""
function initial_conditions!(N::Integer, phi::AbstractVector)
    phi .= 2π *rand(N)
    z = sum(exp.(1im .* phi)) / N

    return abs(z), angle(z)
end

"""
    step!(N::Integer, w::Float64, q::Float64, s::Float64, r::Float64, psi::Float64, dt::Float64, 
    sqdt::Float64, phi::AbstractVector, zk::AbstractVector, nharms::Integer)

Advances by a single dt the simulation. Needs the simulation parameters `N,w,q,s`, the mean-field Kuramoto order parameter at the current 
time `r,psi`. Then, timesteps, and state of the oscillators `phi`. This function also updates `zk`, the value of the Kuramoto-Daido parameters, in place.
`nharms` is the number of KD parameters to compute. The more parameters, the slower it is.

    Return: R, psi
"""
function step!(N::Integer, w::Float64, q::Float64, s::Float64, r::Float64, psi::Float64, dt::Float64, sqdt::Float64, 
    phi::AbstractVector, zk::AbstractVector, nharms::Integer)

    zk .= 0. + 0im 
    @inbounds @simd for i=1:N
        phi[i] += dt * (w + q*r*sin(psi - phi[i])) + sqdt*s*randn()
        phi[i] = mod(phi[i], 2π)

        auxc = cos(phi[i])
        auxs = sin(phi[i])

        oldc = auxc
        olds = auxs

        zk[1] += (oldc + 1im* olds) 
        @simd for k=2:nharms
            nextc = auxc * oldc - auxs * olds 
            nexts = auxs * oldc + auxc * olds 
            zk[k] += (oldc + 1im* olds) 
            oldc = nextc 
            olds = nexts
        end
    end
    
    zk ./= N

    return abs(zk[1]), angle(zk[1])
end


"""
    step_relax(N::Integer, w::Float64, q::Float64, s::Float64, r::Float64, psi::Float64, dt::Float64, sqdt::Float64, phi::AbstractVector)

As step, but it does only update the system state `phi` and returns the Kuramoto-Daido parameter to be able to integrate the next step. 
    Return: R, psi
"""

function step_relax!(N::Integer, w::Float64, q::Float64, s::Float64, r::Float64, psi::Float64, dt::Float64, sqdt::Float64, 
    phi::AbstractVector)
    auxc, auxs = 0., 0.
    @inbounds @simd for i=1:N
        phi[i] += dt * (w + q*r*sin(psi - phi[i])) + sqdt*s*randn()
        phi[i] = mod(phi[i], 2π)

        auxc += cos(phi[i])
        auxs += sin(phi[i])
    end
    z = (auxc + 1im * auxs)/N
    return abs(z), angle(z)
end


"""
    get_phase_diagram(N::Integer, w::Float64, q0::Float64, qf::Float64, nq::Integer, s2::Float64, trelax::Float64, tf::Float64, filename::String; dt=0.01, nsample=100)

Make a phase diagram with a system of `N` oscillators, including frequency `w`, sweeping the coupling from `q0` to `qf`, getting `nq` points,
and having noise instensity squared `s2``.
The system will relax for a time `trelax` and then write the values of the first Kuramoto order and its variance at the file
`filename` for a time `tf`. One can set timestep `dt` and the number of iterations between writing to the file `nsample`.

Output file format: q avg_r(q) var_r(q) [number of rows: nq] 
"""
function get_phase_diagram(N::Integer, w::Float64, q0::Float64, qf::Float64, nq::Integer, s2::Float64, trelax::Float64, tf::Float64, filename::String; 
    dt=0.01, nsample=100)

    #Get constants 
    sqdt = sqrt(dt)
    s = sqrt(s2)

    #Times to integrate. We will split the relaxation time in two.
    nrelax_initial = floor(0.75*trelax/dt)
    nrelax_thermal = floor(0.25*trelax/dt)
    nsteps = floor(tf/dt)

    #Initialize the system
    phi = Vector{Float64}(undef, N)
    r, psi = initial_conditions!(N, phi)

    #Discretization of space
    dq = (qf - q0) / nq

    #Let the system relax at the beginning, starting on the smallest q
    for i=1:nrelax_initial
        r, psi = step_relax!(N, w, q0, s, r, psi, dt, sqdt, phi)
    end

    #Get the phase diagram and put it into a file
    open(filename, "w") do output
        for q=q0:dq:qf

            #Thermalize. Initial state of the simulation is last state of the last ones
            #to minimize required thermalization time
            for i=1:nrelax_thermal
                r, psi = step_relax!(N, w, q, s, r, psi, dt, sqdt, phi)
            end

            #Iniatiate measurements
            avr = 0.
            avr2 = 0.
            nmeasures = 0
                
            #For phase diagram we just want R, so we can always use the step_relax to save some CPU cycles
            for i=1:nsteps
                r, psi = step_relax!(N, w, q, s, r, psi, dt, sqdt, phi)
                
                #Compute the average every nsample iterations
                if i%nsample==1
                    avr += r
                    avr2 += r*r
                    nmeasures += 1
                end
            end

            #Finish the average and write to the file
            avr  /= nmeasures
            avr2 /= nmeasures

            write(output, "$q $avr $(avr2 - avr*avr)\n")
        end
    end


end



"""
    get_timeseries(N::Integer, w::Float64, q::Float64, s2::Float64, trelax::Float64, tf::Float64, filename::String; dt=0.01, nsample=100, nharms=5)

Make a timeseries with a system of `N` oscillators, including frequency `w`, coupling `q`, noise instensity squared `s2``.
The system will relax for a time `trelax` and then write the values of all the Kuramoto-Daido amplitudes to 
`filename` for a time `tf`. One can set timestep `dt` and the number of iterations between writing to the file `nsample`, as well as 
the amount of Kuramoto-Daido parameters to be computed, `nharms`.

Output file format: time r_1(t) r_2(t) ... r_nharm(t) [number of rows: tf/nsample]
"""
function get_timeseries(N::Integer, w::Float64, q::Float64, s2::Float64, trelax::Float64, tf::Float64, filename::String; dt=0.01, nsample=100, nharms=5)
    #Constants
    sqdt = sqrt(dt)
    s = sqrt(s2)

    #Iterations to integrate
    nrelax = floor(trelax/dt)
    nsteps = floor(tf/dt)

    #Initialize system state and get first Kuramoto parameter
    phi = Vector{Float64}(undef, N)
    r, psi = initial_conditions!(N, phi)

    #Relax the system
    for i=1:nrelax
        r, psi = step_relax!(N, w, q, s, r, psi, dt, sqdt, phi)
    end

    #Initialize all the KD
    zk = zeros(nharms) + 1im * zeros(nharms)

    #Open file and start writing simulation results
    open(filename, "w") do output 
        for i=1:nsteps
            r, psi = step!(N, w, q, s, r, psi, dt, sqdt, phi, zk, nharms)

            #EAch nsample iterations, we write the results to a file
            #We write all Kuramoto-Daido amplitudes by columns
            if (i % nsample == 1)
                write(output, "$(i*dt) $r ")
                for k=2:nharms-1
                    write(output, "$(abs(zk[k])) ")
                end
                write(output, "$(abs(zk[nharms]))\n")
            end
        end
    end

    return nothing
end