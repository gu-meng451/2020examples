using Pkg
Pkg.activate(".")
using Plots, LinearAlgebra, Printf

## Problem setup
# Inspired by <https://www.mathworks.com/company/newsletters/articles/stiff-differential-equations.html>
delta = 0.01
f(y,t) = [ y[1]^2 * (1-y[1]) ]
y0 = delta
tf = 2/delta

##
include("BDF.jl")
h = 1
x,t = BDF.bdf4(f, tf, h, y0, tol=1e-8)

plot(t, x[:,1], label=@sprintf("bdf4, h = %.3e", h),
    xlabel="Time t",
    ylabel="Output y(t)" )

## Variable BDF:
# include("vBDF.jl")
# x,t = vBDF.vbdf(f, tf, y0, h, solverTol=1e-8, subStepTol=1e-4)
# 
# plot!(t, x[:,1], label="vbdf3" )

## TODO: add adaptive RK to compare
