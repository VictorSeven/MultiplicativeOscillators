"""
Contains code to integrate the full equations of the Kuramoto-Daido parameters
model using Cartesian coordinates
"""
module KuramotoCartesian
export get_timeseries, phase_diagram

using LinearAlgebra, Distributions, Random

"""
    get_corrs!(nharm::Int, a::Float64, x::Vector, y::Vector, corr::Matrix)

Computes the multiplicative noise matrix and stores it in `corr`, assuming
`nharm` harmonics, noise intensity `a` and Kuramoto-Daido parameters represented
in Cartesian coordinates `x` and `y`.
"""
function get_corrs!(nharm, a, x, y, corr)
    #Fill half of the matrix (a corr matrix is symmetric)
    for i=1:nharm
        for j=1:i-1
            #Closure condition. All harmonics larger than nharm are 0
            if i+j <= nharm 
                xterm = x[i+j] 
                yterm = y[i+j]
            else
                xterm = 0.0 
                yterm = 0.0
            end

            #Matrix element
            corr[j,i] = j*i*a*(x[i-j] - xterm)
            corr[j+nharm,i] = j*i*a*(y[i-j] - yterm)
            corr[j+nharm,i+nharm] =  j*i*a*(x[i-j] + xterm)
        end

        #Fill the diagonal
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

    #Enforce symmetry
    corr .= Symmetric(corr, :L)
end

"""
    get_correlated_vars!(corr::Matrix, xi::Vector; thres=1e-10, maxits=10)

Obtain the stochastic noise from the correlation matrix `corr` and store it in the 
vector `xi`. This method will try to avoid ill-defined correlation matrices by capping
eigenvalues smaller than `thres` to `thres` (ensures positive definite) and symmetrizing 
after eigenvalues are capped. The check is repeated after symmetrization, for a total of
`maxits` iterations.
"""
function get_correlated_vars!(corr, xi; thres=1e-10, maxits=10)
    #Is our matrix positive definite? 
    #Check and overwrite corr with Cholesky decomposition
    good_corrs = isposdef!(corr)
    nits = 0
    #If not,
    while nits < maxits && !good_corrs 
        #Compute the eigenvalues 
        eigsys = eigen(corr)
        #Cap eigenvalues to the threshold and rewrite the matrix in our basis
        corr .= eigsys.vectors * (max.(eigsys.values, thres) .* eigsys.vectors')
        #Symmetrize
        corr .= 0.5 * (corr + corr')
        #Repeat the check (symmetrization might include problems)
        eigsys = eigen(corr)
        good_corrs = !any(eigsys.values .< thres)
        nits += 1
    end

    #Now corr contains the Cholesky matrix, so use it to obtain the random numbers
    randn!(xi)
    xi .= corr * xi
    return nothing
end


"""
    step!(nharm::Int, old_x::Vector, old_y::Vector, x::Vector, y::Vector, corr:Matrix, w::Float64, q::Float64, s2::Float64, a::Float64, dt::Float64, sqdt::Float64, t::Float64, xi::Vector)

An integration step using `nharm` harmonics, computing the new Kuramoto-Daido parameters `x`, `y` from
the old ones `old_x` and `old_y` (cartesian coordinates). Uses the frequency `w`, coupling `q`, noise intensity `s2` 
square and mesoscopic noise `a`. The variables `corr` and `xi` are the correlation matrix and stochastic
noise and are given to avoid allocating memory. The method also needs the timestep `dt`, its 
square root `sqdt` and current timestep `t`.
"""
function step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
    #Get the correlation matrix and generate multiplicative noise
    get_corrs!(nharm, a, old_x, old_y, corr)
    get_correlated_vars!(corr, xi)

    #Computation of the first harmonic 
    #We include 0th harmonic manually, since the 0th is always r[0]=1)
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

"""
    initical_conditions!(nharm::Int, x::Vector, y::Vector)

Set the initial conditions, filling the vectors `x`, `y` (of size `nharm`) of Kuramoto-Daido amplitudes.
"""
function initial_conditions!(nharm, x, y)
    r = Vector{Float64}(undef, nharm)
    psi = Vector{Float64}(undef, nharm)

    r[1] = 0.01*rand()
    psi[1] = 2*π*rand() 
    k = collect(2:nharm)
    r[2:nharm] = r[1] .^ k 
    psi[2:nharm] = psi[1] .* k

    #Get them in cartesian coordinates
    @. x = r * cos(psi) 
    @. y = r * sin(psi) 
end

"""
    get_timeseries(nharm::Int, t_thermal::Float64, tf::Float64, w::Float64, q::Float64, sys_size::Int, s2::Float64, fpath::String; dt=0.01, nsample=10 )

Make a timeseries with a system of `nharm` harmonics, including frequency `w`, coupling `q`, system size `sys_size`
and noise intensity squared `s2`. The system will relax for a time `t_thermal` and then write the 
values of all the Kuramoto-Daido amplitudes to `fpath` for a time `tf`. One can set timestep `dt` and
the number of iterations between writing to the file `nsample`.

Output file format: time r[1] r[2] ... r[nharm]
Observe that the output are the Kuramoto-Daido amplitudes
"""
function get_timeseries(nharm, t_thermal, tf, q, sys_size, s2, fpath; w=0.01, dt=0.01, nsample=10)
    #Get constants we need for the integration
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 

    #Initialize vectors
    old_x = Vector{Float64}(undef, nharm)
    old_y = Vector{Float64}(undef, nharm)
    x = Vector{Float64}(undef, nharm)
    y = Vector{Float64}(undef, nharm)

    corr = zeros(2*nharm, 2*nharm)

    #Initial conditions
    initial_conditions!(nharm, old_x, old_y)

    #Thermalize
    for t=0:dt:t_thermal 
        step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
        old_x, old_y, x, y = x, y, old_x, old_y
    end

    #Start writing the timeseries
    open(fpath, "w") do output 
        t = 0
        nt = 0
        while t < tf 
            step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
            old_x, old_y, x, y = x, y, old_x, old_y

            #Do not write all the iterations to avoid too dense file if dt is small
            #Compute amplitudes R from cartesian
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
function phase_diagram(nharm, t_thermal, tf, q0, qf, nq, sys_size, s2, fpath; sampling=100, w=0.1, dt=0.01)
    #Constants of the system
    a = 0.5 * s2 / sys_size
    sqdt = sqrt(dt) 

    #Initialize the vectors
    q_values = LinRange(q0, qf, nq)

    old_x = Vector{Float64}(undef, nharm)
    old_y = Vector{Float64}(undef, nharm)
    x = Vector{Float64}(undef, nharm)
    y = Vector{Float64}(undef, nharm)

    xi = zeros(2*nharm)
    corr = zeros(2*nharm, 2*nharm)

    open(fpath, "w") do output 
        for q in q_values
            #Initial conditions
            initial_conditions!(nharm, old_x, old_y)

            #Thermalize
            for t=0:dt:t_thermal 
                step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)
                old_x, old_y, x, y = x, y, old_x, old_y
            end

            #Prepare to measure
            avr = 0.0
            avr2 = 0.0
            nmeasures = 0

            t = 0
            nt = 0
            #Measure
            while t < tf 
                step!(nharm, old_x, old_y, x, y, corr, w, q, s2, a, dt, sqdt, t, xi)

                #Here we take a measure, every sampling iterations
                if (nt % sampling == 0)
                    r = sqrt(x[1]^2 + y[1]^2) 
                    avr += r 
                    avr2 += r*r
                    nmeasures += 1
                end

                #Prepare for next iteration
                old_x, old_y, x, y = x, y, old_x, old_y

                t += dt
                nt += 1
            end

            #Finish the average
            avr /= nmeasures
            avr2 /= nmeasures

            write(output, "$q $avr $(avr2-avr*avr)\n")
        end
    end

    return nothing 
end


end #Module end