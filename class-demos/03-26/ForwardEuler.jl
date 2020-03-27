module ForwardEuler

using LinearAlgebra
export feuler

function feuler(f, tf, dt, x0)
    time = 0:dt:tf
    n = length(time)
    p = length(x0)
    x = zeros(n,p)

    x[1,:] .= x0

    for i = 1:n-1
        x[i+1,:] = x[i,:] + dt*f(x[i,:], time[i])
    end

    return x, time

end

end
