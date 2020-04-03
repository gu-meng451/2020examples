using Pkg
Pkg.activate(".")

using Plots

## Load the Backwards Euler Method:
include("../3-26-20/BackwardEuler.jl")
using .BackwardEuler

## Define the Problem
λ = -2.0
f(x,t) = λ*x
tf = 4.0
x0 = 1.0
h  = 0.1

x_exact(t) = exp.(λ*t)*x0

plot( x_exact, 0, tf, label="Exact Solution",
  xlabel="Time t [s]", ylabel="x(t)")

## Backward Euler
x,time = beuler(f, tf, h, x0)
plot!( time, x, label="Backward Euler" )

## Compare to Foward Euler
include("../03-26/ForwardEuler.jl")
using .ForwardEuler

x,time = feuler( f, tf, h, x0)

plot!( time, x, label="Forward Euler" )

##
include("AB2.jl")
using .AB2

x,time = ab2(f, tf, h, x0)
plot!(time, x, label="AB2")

##
include("AM3.jl")
using .AM3

x,time = am3(f, tf, h, x0)
plot!( time, x, label="AM3")
