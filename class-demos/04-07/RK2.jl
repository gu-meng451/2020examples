module RK2

using LinearAlgebra


function rk2_step(f, tn, h, xn; alpha = 2/3)

    c1 = alpha
    a21 = alpha
    b1 = 1 - 1/2/alpha
    b2 = 1/2/alpha

    k1 = f(xn, tn)
    k2 = f(xn + h*a21*k1, tn + h*c1 )

    xn_p1 = xn + h*(b1*k1 + b2*k2)

    return xn_p1

end

function rk2( f, tf, h, x0 )

    time = 0:h:tf
    n = length(time)
    p = length(x0)

    x = zeros(n,p)

    # Initial Condition
    x[1,:] .= x0

    # loop over the rest of time
    for i = 1:n-1
        x[i+1,:] = rk2_step(f, time[i], h, x[i,:] )
    end

    return x, time
end


function rk2_vstep(f, tn, h, xn; alpha=2/3, tol = 1e-5)

        x1  = rk2_step(f, tn, h, xn, alpha = alpha)

        xs1 = rk2_step(f, tn, h/2, xn, alpha = alpha)
        xs  = rk2_step(f, tn+h/2, h/2, xs1, alpha=alpha)

        return

    end


end
