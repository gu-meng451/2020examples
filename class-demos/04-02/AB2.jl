module AB2

using LinearAlgebra

function ab2(f, tf, h, x0)
    time = 0:h:tf
    n = length(time)
    p = length(x0)
    x = zeros(n,p)

    # Initial Condition
    x[1,:] .= x0

    # 1 Euler Step (this is cheap and possibly bad)
    fn = f(x[1,:], time[1])
    x[2,:] = x[1,:] + h*fn

    # AB2 the rest
    for i = 2:n-1
        fn_m1 = fn # f(x[i-1,:], time[i-1])
        fn = f(x[i,:], time[i])
        #x[i+1,:] = x[i,:] + h/2*( 3*f(x[i,:], time[i]) - f(x[i-1,:], time[i-1])  )
        x[i+1,:] = x[i,:] + h/2*( 3*fn - fn_m1 )
    end

    return x, time
end

end
