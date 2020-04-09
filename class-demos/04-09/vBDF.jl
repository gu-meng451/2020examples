module vBDF

using LinearAlgebra
using ForwardDiff

export vbdf

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


function trap_step(f, xn, tn, h, xstar, tol, iterMax)

    # Residual equation
    g(y) = -y + xn + h/2*( f(y,tn+h) + f(xn,tn) )

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function trap_adaptStep(f, xn, tn, hn, solverTol, iterMax, maxSubSteps, subStepTol)

    h = hn
    substep = 0
    flag = 0

    while flag == 0
        substep += 1

        # Single step to h
        ys = forwardEuler_step(f, xn, tn, h)
        y  = trap_step(f, xn, tn, h, ys, solverTol, iterMax)

        # 2 Steps top h
        z1s = forwardEuler_step(f, xn, tn, h/2)
        zm  = trap_step(f, xn, tn, h/2, z1s, solverTol, iterMax)
        z2s = forwardEuler_step(f, zm, tn+h/2, h/2 )
        z   = trap_step(f, zm, tn+h/2, h/2, z2s, solverTol, iterMax)

        err = norm( z - y )
        if err <= subStepTol
            return z , h, err
        elseif substep >= maxSubSteps
            error("Took too many substeps and failed to converge")
        else
            h = h/2
        end

    end

end

function bdf2_adaptStep(f, xn, xn_m1, tn,  tn_m1, hn, solverTol, iterMax, maxSubSteps, subStepTol)

    h = hn
    substep = 0
    flag = 0

    while flag == 0
        substep += 1

        # Prediction ystar
        xs = ab2_vstep(f, h, xn, xn_m1, tn, tn_m1)

        # O(h^2) corrector
        xn_p1 = bdf2_vstep(f, h, xn, xn_m1, tn, tn_m1, xs, solverTol, iterMax)

        # O(h) Corrector
        y = backwardEuler_step(f, xn, tn, h, xs, solverTol, iterMax)

        err = norm( xn_p1 - y )

        if err <= subStepTol
            return xn_p1, h, err
        elseif substep >= maxSubSteps
            error("Took too many substeps and failed to converge")
        else
            h = h/2
        end

    end

end

function bdf3_adaptStep(f, xn, xn_m1, xn_m2, tn, tn_m1, tn_m2, hn, solverTol, iterMax, maxSubSteps, subStepTol)

    h = hn
    substep = 0
    flag = 0

    while flag == 0
        substep += 1

        # Prediction ystar
        xs = ab2_vstep(f, h, xn, xn_m1, tn, tn_m1)

        # O(h^2) corrector
        y = bdf2_vstep(f, h, xn, xn_m1, tn, tn_m1, xs, solverTol, iterMax)

        # O(h^3) Corrector
        xn_p1 = bdf3_vstep(f, h, xn, xn_m1, xn_m2, tn, tn_m1, tn_m2, y, solverTol, iterMax)

        err = norm( xn_p1 - y )

        if err <= subStepTol
            return xn_p1, h, err
        elseif substep >= maxSubSteps
            error("Took too many substeps and failed to converge")
        else
            h = h/2
        end

    end

end

function ab2_vstep(f, h, xn, xn_m1, tn, tn_m1)
    h_m1 = tn - tn_m1
    bn_m1 = -h^2/2/h_m1
    bn    = h + h^2/2/h_m1

    return xn + bn*f(xn,tn) + bn_m1*f(xn_m1,tn_m1)
end

function bdf2_vstep(f, h, xn, xn_m1, tn, tn_m1, xstar, tol, iterMax)
    h0 = tn - tn_m1
    a1 = h/(h0*(h + h0))
    a2 = -(h + h0)/(h*h0)
    a3 = (2*h + h0)/(h*(h + h0))

    # Residual equation
    g(y) = a3*y + a2*xn + a1*xn_m1 - f(y,tn+h)

    # Solve:
    x1, iter = newtonsolve(g, xstar, tol, iterMax)

    return x1

end

function bdf3_vstep(f, h, xn, xn_m1, xn_m2, tn, tn_m1, tn_m2, xstar, tol, iterMax)

    h1 = tn - tn_m1
    h0 = tn_m1 - tn_m2

    a4 = 1/(h + h0 + h1) + 1/(h + h1) + 1/h
    a3 = -(h^2 + h*h0 + 2*h*h1 + h0*h1 + h1^2)/(h*h1*(h0 + h1))
    a2 = h*(h + h0 + h1)/(h0*h1*(h + h1))
    a1 = -h*(h + h1)/(h0*(h*h0 + h*h1 + h0^2 + 2*h0*h1 + h1^2))

    # Residual equation
    g(y) = a4*y + a3*xn + a2*xn_m1 + a1*xn_m2 - f(y,tn+h)

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

function vbdf(f, tf, x0, h0; solverTol = 1e-4, solverIterMax=4*length(x0), maxSteps = 10*round(Int, tf/h0), subStepTol=1e-6, maxSubSteps = 20 )
    time = []
    p = length(x0)
    x = Array{Float64,2}(undef,0,p)

    ## Initial Condition
    push!(time, 0. )
    x = [x; x0']

    h = h0
    j = 1

    # TODO: count fails
    fail_count = 0

    for i = 1:maxSteps

        if i == 1
            y,h,err = trap_adaptStep(f, x[i,:], time[i], h, solverTol, solverIterMax, maxSubSteps, subStepTol)

        elseif i == 2
            y,h,err = bdf2_adaptStep(f, x[i,:], x[i-1,:], time[i], time[i-1], h, solverTol, solverIterMax, maxSubSteps, subStepTol)

        else
            y,h,err = bdf3_adaptStep(f, x[i,:], x[i-1,:], x[i-2,:], time[i], time[i-1], time[i-2], h, solverTol, solverIterMax, maxSubSteps, subStepTol)

        end

        # Add the new point in time
        x = vcat( x, y' )
        push!(time, time[i]+h )

        # predict a good next time step
        h = 0.9*h*min(max(1/err, 0.3),2)

        if time[end] >= tf
            break
        end

    end

    return x, time

end


end
