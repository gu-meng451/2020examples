# Julia version of Matlab code:
using Pkg
Pkg.activate(".")
using Plots, Printf
using DifferentialEquations

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

# solve the problem
sol = solve(prob)

# Let's make an example plot:
plt1 = plot(sol,
        linewidth=2,
        xaxis="Time t [s]",
        # label="\\theta_1 [rad]",
        yaxis="\\theta_1 [rad]",
        label="",
        vars = (0,1) );
plt2 = plot(sol,
        linewidth=2,
        xaxis="Time t [s]",
        yaxis="\\theta_2 [rad]",
        label="",
        vars = (0,2) );
plot(plt1, plt2, layout=(2,1) )

# We can also build an animation.  I'll make some functions (just like matlab)
# that handle the interpolations at any time t.  It's basically `deval` built
# right in.
θ1(t) = sol(t)[1]
θ2(t) = sol(t)[2]
x1(t) = l₁*sin(θ1(t))
y1(t) = -l₁*cos(θ1(t))
x2(t) = x1(t) + l₂*sin(θ2(t))
y2(t) = y1(t) - l₂*cos(θ2(t))

nframes = 50*tf
time = LinRange(0, tf, round(Int,nframes) )

# I'm going to use some fancier stuff, so I'll switch the backend of Plots.jl
# to be PyPlot for the animation.
pyplot()
plt = plot([0,x1(0)], [0,y1(0)], lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(-l₁-l₂, +l₁+l₂),
        ylim=(-l₁-l₂, +l₁+l₂),
        dpi=300 );
plot!(plt, [x1(0), x2(0)], [y1(0), y2(0)], lw=3, label="");
plot!(plt, title=@sprintf("Time t = %.3f [s]", 0),
        titlefont=font("DejaVu Sans Mono") )
anim = @animate for t in time
    plt[1] = ( [0,x1(t)], [0,y1(t)] )
    plt[2] = ([x1(t), x2(t)], [y1(t), y2(t)] )
    plot!(plt, title=@sprintf("Time t = %6.3f [s]", t) )
end
gif(anim, "test1.gif", fps=20)
