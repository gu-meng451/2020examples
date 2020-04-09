using Pkg
Pkg.activate(".")
using Plots, LinearAlgebra, Printf

ζ = 0.3
x0 = [1, 1]
tf = 10.

f(z,t) = [ z[2], -z[1] - 2*ζ*z[2] ]

## Exact solution
x_exact(t) = exp(-ζ*t)*( x0[1]*cos( sqrt(1-ζ^2)*t ) + (x0[2] + ζ*x0[1])/sqrt(1-ζ^2)*sin(sqrt(1-ζ^2)*t) )

##
h0 = 0.1
include("vBDF.jl")
x,t = vBDF.vbdf(f, tf, x0, h0, subStepTol=1e-3)

plt = plot(x_exact, 0, tf, label="exact",
    lw = 3,
    xlabel="Time t", ylabel="output x(t)")
plot!(plt, t, x[:,1], linestyle=:dash, marker=:circle, label="vBDF")
