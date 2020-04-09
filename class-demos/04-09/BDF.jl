module BDF

using LinearAlgebra
using ForwardDiff


function forwardEuler_step(f, xn, tn, h)
    return xn + h*f(xn, tn)
end

function backwardEuler_step(f, xn, tn, h, xstar, tol, iterMax)

    # Residual equation
    g(y) = -y + xn + h*f(y,tn+h)

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function ab2_step(f, xn, xn_m1, tn, h)
    return xn + h/2*( 3*f(xn,tn) - f(xn_m1,tn-h) )
end

function bdf2_step(f, xn, xn_m1, tn, h, xstar, tol, iterMax)

    # Residual equation
    g(y) = y - 4/3*xn + 1/3*xn_m1 - 2/3*h*f(y,tn+h)

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function bdf3_step(f, xn, xn_m1, xn_m2, tn, h, xstar, tol, iterMax)

    # Residual equation
    g(y) = y - 18/11*xn + 9/11*xn_m1 - 2/11*xn_m2 - 6/11*h*f(y,tn+h)

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function bdf4_step(f, xn, xn_m1, xn_m2, xn_m3, tn, h, xstar, tol, iterMax)

    # Residual equation
    g(y) = y - 48/25*xn + 36/25*xn_m1 - 16/25*xn_m2 + 3/25*xn_m3 - 12/25*h*f(y,tn+h)

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function newtonsolve(g, ys, tol, iterMax)
    # compute derivative
    ∇g(y) = ForwardDiff.jacobian( g, y)

    # Newton's method
    flag = 0
    iter = 0
    y0 = copy(ys)
    g0 = g(y0)
    while flag == 0
        iter += 1
        y1 = y0 - ∇g(y0)\g0
        y0[:] = y1
        g0 = g(y0)
        if norm( g0 ) <= tol
            flag = 1
        elseif iter >= iterMax
            error("Newton's Method Failed to Converge")
        end
    end
    return y0, iter
end

function bdf4(f, tf, h, x0; tol=1e-6, iterMax=4*length(x0) )
    time = 0:h:tf
    n = length(time)
    p = length(x0)
    x = zeros(n,p)

    ## Initial Condition
    x[1,:] .= x0

    for i = 1:n-1
        if i == 1
            ## Step 1: Euler
            xs = forwardEuler_step(f, x[i,:], time[i], h)
            x[i+1,:] = backwardEuler_step(f, x[i,:], time[i], h, xs, tol, iterMax)
        elseif i == 2
            ## step 2:
            xs = ab2_step(f, x[i,:], x[i-1,:], time[i], h)
            x[i+1,:] = bdf2_step(f, x[i,:], x[i-1,:], time[i], h, xs, tol, iterMax)

        elseif i == 3
            ## step 3:
            xs = ab2_step(f, x[i,:], x[i-1,:], time[i], h)
            x[i+1,:] = bdf3_step(f, x[i,:], x[i-1,:], x[i-2,:], time[i], h, xs, tol, iterMax)
        else
            xs = ab2_step(f, x[i,:], x[i-1,:], time[i], h)
            x[i+1,:] = bdf4_step(f, x[i,:], x[i-1,:], x[i-2,:], x[i-3,:], time[i], h, xs, tol, iterMax)
        end

    end

    return x, time

end



end
