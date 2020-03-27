using Pkg
Pkg.activate(".")

using Plots
include("ForwardEuler.jl")
using .ForwardEuler

## Define the function to integrate
λ = -2
f(x,t) = λ*x

x0 = 1.0
h  = 0.1
tf = 4.0

## Solve the system
X,T = feuler(f, tf, h, x0)

##
x_exact(t) = exp.(λ*t)*x0

plot( x_exact, 0, tf, label="Exact Solution" )
plot!( T, X, label="ForwardEuler")
plot!(xlabel="Time t [s]",
    ylabel="x(t)")
