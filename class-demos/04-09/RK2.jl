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

        return xs, x1

end


function rk2v( f, tf, h0, x0; Atol = 1e-4, Rtol = 1e-4, maxSteps = 10*round(Int, tf/h0) )

    time = [0.]
    x = x0'
    j = 1
    h = h0
    numberOfFailedSteps = 0

    for i = 1:maxSteps

        # xs, x1 = rk2_vstep(f, time[j], h, x[j,:])
        xs, x1 = rk_BS(f, time[j], h, x[j,:])

        # compute error:
        err = compute_err(Atol, Rtol, x[j,:], x1, xs)

        if err <= 1
            # step was good
            j = j+1
            push!(time, time[end]+h)
            x = vcat(x, xs')
        else
            numberOfFailedSteps += 1
        end

        if time[end] >= tf
            break
        end

        h = h*min( 2., max(0.3, (1/err)^(1/(1+2)) ))

    end

    return x, time, numberOfFailedSteps


end

function compute_err(Atol, Rtol, xn, x1, xs)
    n = length(xn)

    sc = [ Atol .+ max(abs(xn[i]), abs(x1[i]))*Rtol for i = 1:n ]

    err = sqrt(  sum( ((x1-xs)./sc).^2 )/n  )

    return err

end


function rk_BS(f, tn, h, xn )
    #Bogacki-Shampine Method 3-2
    c1 = 0.
    c2 = 1/2
    c3 = 3/4
    c4 = 1
    a21 = 1/2
    a31 = 0. # hi
    a32 = 3/4
    a41 =2/9
    a42 =1/3
    a43 =4/9
    b1 = 2/9
    b2 = 1/3
    b3 = 4/9
    b4 = 0
    bt1 = 7/24
    bt2 = 1/4
    bt3 = 1/3
    bt4 = 1/8

    k1 = f(xn, tn+c1*h)
    k2 = f( xn + h*a21*k1, tn+c2*h )
    k3 = f( xn + h*a31*k1 + h*a32*k2, tn+c3*h )
    k4 = f( xn + h*a41*k1 + h*a42*k2 + h*a43*k3, tn+c4*h)

    # O(h^3)
    xs = xn + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4)

    # O(h^2)
    x1 = xn + h*(bt1*k1 + bt2*k2 + bt3*k3 + bt4*k4)

    return xs, x1
end



end
