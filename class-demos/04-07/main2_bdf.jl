using Pkg
Pkg.activate(".")
using Plots, LinearAlgebra, Printf

ζ = 0.3
x0 = [0.1, 0.2]
tf = 10.

## Exact solution
x_exact(t) = exp(-ζ*t)*( x0[1]*cos( sqrt(1-ζ^2)*t ) + (x0[2] + ζ*x0[1])/sqrt(1-ζ^2)*sin(sqrt(1-ζ^2)*t) )
f(z,t) = [z[2], -z[1] - 2*ζ*z[2]  ]

##
include("BDF.jl")
h = 0.1
xs = BDF.forwardEuler_step(f, x0, 0, h)
x1 = BDF.backwardEuler_step(f,x0,0,h,xs, 1e-6, 200)

x2 = BDF.bdf2_step(f, x1, x0, 0+h, h, x1, 1e-6, 200)
