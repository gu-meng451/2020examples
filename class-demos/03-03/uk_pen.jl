using Pkg
Pkg.activate(".")
using DifferentialEquations, Plots,
      LinearAlgebra, ForwardDiff, Printf

## Setup our parameters
m = 1.
l = 0.5
g = 9.81

tf = 10.

# Initial conditions in terms of the angle theta:
θ0 = 45*pi/180
θd0 = 0.

## Let's build the Udwadia-Kalaba Equation
# q = [x, y]

# Mass Matrix
M = diagm(0=>[m,m])

# External Forcing
Qex = [0, -m*g]

# The constraint:
ϕ(q) = q[1]^2 + q[2]^2 - l^2

# Build the U-K Equation
A(q) = ForwardDiff.gradient(ϕ, q)'
b(q,qd) = - qd'*ForwardDiff.hessian(ϕ,q)*qd
Ms = sqrt(M) |> real

qdd(q,qd)= M\( Qex + Ms*pinv(A(q)/Ms)*( b(q,qd) - (A(q)/M)*Qex) )

## Setup and Solve the ODE IVP:
function myode!( dz, z, p, t)
    q = z[1:2]
    qd= z[3:4]

    dz[1:2] = qd
    dz[3:4] = qdd(q,qd)
end

# Initial Conditions:
z0 = [l*sin(θ0), -l*cos(θ0), l*θd0*cos(θ0), l*θd0*sin(θ0)]

# define the 'prob'
tspan = (0., tf)
p = []
prob = ODEProblem(myode!,z0,tspan,p)

# solve (using defaults)
sol = solve(prob)

## Solve
plot(sol)

## Plot the length
plt_L = plot( t->ϕ(sol(t)), 0, tf )


## Let's make an animation
x(t) = sol(t)[1]
y(t) = sol(t)[2]

t1 = LinRange(0,tf,500)
xmin = minimum(x,t1)
xmax = maximum(x,t1)
ymin = minimum(y,t1)
ymax = max(0, maximum(y,t1) )


nframes = 50*tf
time = LinRange(0, tf, round(Int,nframes) )
plt = plot([0,x(0)], [0,y(0)], lw=3, label="",
        aspect_ratio=:equal,
        xaxis="x [m]", yaxis="y [m]",
        xlim=(xmin, xmax),
        ylim=(ymin, ymax),
        dpi=300 )
plot!(plt, title=@sprintf("Time t = %.3f [s]", 0),
        titlefont=font("DejaVu Sans Mono") )

anim = @animate for t in time
    plt[1] = ( [0,x(t)], [0,y(t)] )
    plot!(plt, title=@sprintf("Time t = %6.3f [s]", t) )
end
gif(anim, "test1.gif", fps=20)
