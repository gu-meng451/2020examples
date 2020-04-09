using Pkg
Pkg.activate(".")
using Plots, LinearAlgebra, Printf

ζ = 0.3
x0 = [0.1, 0.2]
tf = 10.

f(z,t) = [ z[2], -z[1] - 2*ζ*z[2] ]

## Exact solution
x_exact(t) = exp(-ζ*t)*( x0[1]*cos( sqrt(1-ζ^2)*t ) + (x0[2] + ζ*x0[1])/sqrt(1-ζ^2)*sin(sqrt(1-ζ^2)*t) )

##
include("BDF.jl")
h = 0.1

x,t = BDF.bdf4(f, tf, h, x0)

plot(x_exact, 0, tf, label="exact")

plot!(t, x[:,1], label="bdf4")
