module TheoryFormulas

using Statistics

export r_6th, r_2th_cumulant, r_oa, finite_size_r, integrate_tyulkina, integrate_full

function r_6th(q, s2)
    ratio = @. s2 * sqrt(2) / q
    root = @. sqrt(complex(51*q*q - 132 * q * s2 + 306*s2*s2)) 
    numerator = @. complex(24*s2 - 9*q - root)
    denominator = @. complex(q - 9*s2)
    return @. real(ratio * sqrt(numerator / denominator))
end

function r_2th_cumulant(q, s2)
    root = @. sqrt(complex(4*q*q + 4*q*s2 - 7*s2*s2))
    numerator = @. complex(2*q - 3*s2 + root) 
    return @. real(0.5 * sqrt(numerator / q))
end

function r_oa(q, s2)
    return @. real(sqrt(complex((q-s2)/q))) 
end

function finite_size_r(q, s2, n)
    aux = @. q * n - s2 * (n - 1)  
    aux = @. aux + sqrt(4*n*q*s2 + aux^2)
    return @. sqrt(aux / (2*q*n))
end

function tyulkina(w, q, s2, oldz, oldx, dt)
    z = oldz + dt * ( 0.5*(q-s2+2*im*w)*oldz - 0.5*q*(oldz * abs2(oldz) + conj(oldz)*oldx)) 
    x = oldx + dt * ( 2*im*w*oldx - s2*(2*oldx + oldz^2) - 2*q*oldx*abs2(oldz) )
    return z,x
end

function full_system(nharm, w, q, s2, oldz, z, dt)
    k = 1 #z[0]=1
    z[k] = oldz[k] + dt *( oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*(oldz[1] - conj(oldz[1])*oldz[k+1]))
    @simd for k=2:nharm-1
        z[k] = oldz[k] + dt *( oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*(oldz[1]*oldz[k-1] - conj(oldz[1])*oldz[k+1]))
    end
    k = nharm #z[nharm+1]=0
    z[k] = oldz[k ] + dt * (oldz[k] * (im*w*k - 0.5*k*k*s2) + 0.5*q*k*oldz[1]*oldz[k-1] )
end

function integrate_tyulkina(w, q, s2; dt=0.01, tf=1000.0, ntrials=10)
    avr = 0.0 
    @fastmath @inbounds for i=1:ntrials
        angles = 2π*rand(100000)
        oldz = mean(exp.(im*angles))
        z2 = mean(exp.(2*im*angles))
        oldx = z2 - oldz.^2  

        z = 0.0
        for t=0.0:dt:tf 
            z, x = tyulkina(w, q, s2, oldz, oldx, dt)
            oldz, oldx = z, x
        end
        avr += abs(z)
    end

    return avr/ntrials 
end

function integrate_tyulkina_simple(w, q, s2; dt=0.01, n=100000, tol=1e-4)
    angles = 2π*rand(n)
    oldz = mean(exp.(im*angles))
    z2 = mean(exp.(2*im*angles))
    oldx = z2 - oldz.^2  
    rnew = abs(oldz)

    eps = 10*tol
    z = 0.0+0.0im 
    for i=1:500000
    #while eps > tol
        z, x = tyulkina(w, q, s2, oldz, oldx, dt)
        #rold = rnew
        #rnew = abs(z)
        #eps = abs((rnew-rold)/rnew)
        oldz, oldx = z, x
    #end
    end

    return abs(oldz)
end

function diagram_tyulkina(q0, qf, nq, w, s2, path; tol=1e-4)
    open(path, "w") do output
        for q in LinRange(q0, qf, nq)
            r = integrate_tyulkina_simple(w, q, s2; tol=tol)
            write(output, "$q $r\n")
        end
    end
end

function timeseries_tyulkina(w, q, s2, outputpath; dt=0.01, tf=1000.0)
    angles = 2π*rand(100000)
    oldz = mean(exp.(im*angles))
    z2 = mean(exp.(2*im*angles))
    oldx = z2 - oldz.^2  

    z = 0.0
    nits = 0
    open(outputpath, "w") do output
        for t=0.0:dt:tf 
            z, x = tyulkina(w, q, s2, oldz, oldx, dt)
            oldz, oldx = z, x
            if nits%10==0
                write(output, "$t $(abs(oldz))\n")
            end
            nits += 1
        end
    end
end

function integrate_full(w, q, s2; dt=0.01, tf=1000.0, ntrials=1, nharm=20)
    avr = 0.0 
    @fastmath @inbounds for i=1:ntrials
        angles = 2π*rand(100000)
        oldz = [mean(exp.(im*angles*k)) for k=1:nharm]
        z = Vector{ComplexF64}(undef, nharm)
        for t=0.0:dt:tf 
            full_system(nharm, w, q, s2, oldz, z, dt)
            oldz, z = z, oldz 
        end
        avr += abs(oldz[1])
    end

    return avr/ntrials 
end

function check_angles(w, q, s2; dt=0.01, tf=1000.0, nharm=30)
    angles = 2π*rand(100000)
    oldz = [mean(exp.(im*angles*k)) for k=1:nharm]
    z = Vector{ComplexF64}(undef, nharm)
    for t=0.0:dt:tf 
        full_system(nharm, w, q, s2, oldz, z, dt)
        oldz, z = z, oldz 
    end
    return angle.(oldz) 
end

end