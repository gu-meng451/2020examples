# using Pkg
# Pkg.activate(".")
using Plots
include("BackwardEuler.jl")
# using .BackwardEuler

λ = -2.0
f(x,t) = λ*x
tf = 4.0
x0 = 1.0
h  = 0.1

##
x,time = BackwardEuler.beuler(f, tf, h, x0)

plot( time, x )
