# Testing in Julia

using Pkg
Pkg.activate(".")
# Pkg.add("BenchmarkTools")
using BenchmarkTools

include("Exp0.jl")
using .Exp0 # this exports Exp into the workspace

x = 10*(2*rand()-1)

@benchmark Exp(x)

@benchmark exp(x)
