module BackwardEuler

using LinearAlgebra

export beuler

function beuler(f, tf, h, x0)

    # x(i+1) = x(i) + h*f(i+1)

    time = 0:h:tf
    n    = length(time)
    p    = length(x0)

    x = zeros(n,p)
    x[1,:] .= x0

    for i = 1:n-1
        # x(i+1) = x(i) + h*f( x(i+1) )
        # y = g(y)
        # Two options: Newton's Method if we know g'(y), fixed point iteration (we'll do this one)

        # Predict via forward Euler x[i+1]
        y = x[i,:] + h*f( x[i,:], time[i] )

        tol = 1e-4
        flag = 0
        iterMax = 200
        iter = 0
        while flag == 0
            iter += 1
            y = x[i,:] + h*f(y, time[i+1])

            residual = norm( y - x[i,:] - h*f(y, time[i+1]) )

            if residual <= tol
                flag = 1
            elseif iter >= iterMax
                flag = -1
                error("Error: failed to converge")
            end
        end

        x[i+1,:] = y

    end

    return x, time

end

end
