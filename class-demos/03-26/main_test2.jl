using Pkg
Pkg.activate(".")

using Plots, BenchmarkTools

## Load the Backwards Euler Method:
include("BackwardEuler.jl")
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

beulerBmark = @benchmark beuler(f, tf, h, x0)

plot!( time, x, label="Backward Euler" )

## Compare to Foward Euler
include("ForwardEuler.jl")
using .ForwardEuler

x,time = feuler( f, tf, h, x0)

feulerBmark = @benchmark feuler(f, tf, h, x0)
# Note how much less memory this takes.  On my laptop it's also 5x faster.

plot!( time, x, label="Forward Euler" )
