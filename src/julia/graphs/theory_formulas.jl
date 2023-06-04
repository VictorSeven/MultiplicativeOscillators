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