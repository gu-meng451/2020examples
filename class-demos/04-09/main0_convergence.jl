using Pkg
Pkg.activate(".")
using Plots, LinearAlgebra, Printf

ζ = 0.3
x0 = [0.1, 0.2]
tf = 10.

## Exact solution
x_exact(t) = exp(-ζ*t)*( x0[1]*cos( sqrt(1-ζ^2)*t ) + (x0[2] + ζ*x0[1])/sqrt(1-ζ^2)*sin(sqrt(1-ζ^2)*t) )
f(z,t) = [z[2], -z[1] - 2*ζ*z[2]  ]

plt1 = plot( t->x_exact(t), 0, tf, label="Exact solution" )
plot!(plt1, xlabel="Time t\\omega_n", ylabel="Output x(t)", grid=true)


## Setup the error function
function err(x,t)
    return norm( x[:,1] - x_exact.(t) )
end

function addXtoPlot!(plt, x, t, name, h)
    plot!(plt, t, x[:,1], label=@sprintf("%s, h = %.3e", name, h) )
end

##
include("../04-02/AM3.jl")
h = 0.1
x,t = AM3.am3(f, tf, h, x0)
addXtoPlot!(plt1, x, t, "AM3", h)

## Let's run this over a range of h:
h_range = 2 .^-(2.0:16.0)

errAM3 = Float64[]
for h in h_range
    x,t = AM3.am3(f, tf, h, x0, tol=1e-8)
    push!(errAM3, err(x,t) )
end
plt2 = plot(h_range, errAM3, marker=:square,
    label="AM3",
    xscale=:log10, yscale=:log10,
    xlabel="Step size h",
    ylabel="Norm of error |x - x_true|_2",
    legend=:bottomright)

## Add AB2
include("../04-02/AB2.jl")
errAB2 = Float64[]
for h in h_range
    x,t = AB2.ab2(f, tf, h, x0)
    push!(errAB2, err(x,t) )
end

plot!(plt2, h_range, errAB2, marker=:star5, label="AB2")

## Forword Euler:
include("../03-26/ForwardEuler.jl")
errFE = Float64[]
for h in h_range
    x,t = ForwardEuler.feuler(f, tf, h, x0)
    push!(errFE, err(x,t) )
end
plot!(plt2, h_range, errFE, marker=:hexagon, label="Forward Euler")

## Backward Euler:
include("../03-26/BackwardEuler.jl")
errBE = Float64[]
for h in h_range
    x,t = BackwardEuler.beuler(f, tf, h, x0, tol=1e-6)
    push!(errBE, err(x,t) )
end
plot!(plt2, h_range, errBE, marker=:star4, label="Backward Euler")


## RK2
include("RK2.jl")
h = 0.1
# x,t = RK2.rk2(f, tf, h, x0)
# addXtoPlot!(plt1, x, t, "RK2", h)

errRK2 = Float64[]
for h in h_range
    x,t = RK2.rk2(f, tf, h, x0)
    push!(errRK2, err(x,t) )
end
plot!(plt2, h_range, errRK2, marker=:star3, label="RK-2")


##
h0 = 0.1
include("RK2.jl")
x,t = RK2.rk2v( f, tf, h0, x0)

plot!(plt1, t, x[:,1], marker=:circle, xlabel="Time t", ylabel="output x(t)" )
