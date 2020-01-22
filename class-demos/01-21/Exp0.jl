module Exp0

export Exp

function Exp(n::Int)
    y = 1.0
    e = 2.7182818284590452
    for i = 1:n
        y *= e
    end
    return y
end

function Exp(x::Float64; tol=1e-12)

    flag_InvExp = false
    if x < 0
        x = -x
        flag_InvExp = true
    end

    # Break apart x into n + z
    n = floor(Int,x)
    z = x - n

    # Iterate to compute e^z
    maxIter = 100
    flag = 0
    i = 0
    expZ = 1.0
    newterm = 1.0
    c = 0.0
    while flag == 0
        i += 1
        newterm *= z/i

        # using Kahan summation:
        # Z += newterm
        y = newterm - c
        t = expZ + y
        c = (t-expZ) - y
        expZ = t

        rel_err = abs(newterm)
        if rel_err <= tol
            flag = 1
        elseif i >= maxIter
            flag = -1
        end
    end

    # compute exp(x) = exp(n)*exp(z)
    # since n::Int, then this calls the Int version of Exp
    y = Exp(n)*expZ

    # flip the output
    if flag_InvExp
        y = 1/y
    end

    # return the output
    #return y,i
    return y

end


end
