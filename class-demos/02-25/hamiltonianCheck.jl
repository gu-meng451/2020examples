# Julia version of Matlab code:
using Pkg
Pkg.activate(".")
using Plots
using DifferentialEquations
pyplot()

# Define the parameters
l₁ = 0.5 #[m]
l₂ = 0.5 #[m]
m₁ = 1.  #[kg]
m₂ = 1. #[kg]
g  = 9.81 #[m/s^2]

tf = 20. #[s]

# load our ODE function
include("doublepen.jl")

# Initial conditions (in terms of angles)
θ₁0  = 45*pi/180 #[rad]
θ₂0  = 60*pi/180 #[rad]
dθ₁0 = 10. # [rad/s]
dθ₂0 = -10. #[rad/s]

# Let's transform these coordinate-based IC to our states
p₁0 = l₁^2*m₁*dθ₁0 + m₂*(2*l₁^2*dθ₁0 + 2*l₁*l₂*cos(θ₁0 - θ₂0)*dθ₂0)/2
p₂0 = m₂*(2*l₁*l₂*cos(θ₁0 - θ₂0)*dθ₁0 + 2*l₂^2*dθ₂0)/2

# Setup the ODE problem
p = l₁, l₂, m₁, m₂, g
tspan = (0., tf)
z0 = [ θ₁0, θ₂0, p₁0, p₂0 ]
prob = ODEProblem( doublepen!, z0, tspan, p )

# solve the problem using the defaults
sol1 = solve(prob)

# Let's see how the Hamiltonian is doing
plt_ham = plot( t -> hamiltonian(sol1,p,t)/hamiltonian(sol1,p,0), 0, tf,
        label="default",
        xaxis="Time t [s]",
        yaxis="Normalized Hamiltonian H(t)/H(0)")

sol2 = solve(prob,RadauIIA5() )
plot!(plt_ham, t -> hamiltonian(sol2,p,t)/hamiltonian(sol2,p,0), 0, tf,
        label="RadauIIA5")

sol3 = solve( prob, TRBDF2() )
plot!(plt_ham, t -> hamiltonian(sol3,p,t)/hamiltonian(sol3,p,0), 0, tf,
        label="TRBDF2")

sol4 = solve( prob, SSPSDIRK2(), dt=0.01 )
plot!(plt_ham, t -> hamiltonian(sol4,p,t)/hamiltonian(sol4,p,0), 0, tf,
        label="SSPSDIRK2")
