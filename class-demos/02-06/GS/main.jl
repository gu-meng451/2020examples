
using LinearAlgebra
using Printf

nx = 5
ny = nx

# Make the array of T, and assume the initial guess is zero.
T = zeros(nx+2, ny+2)


# Apply the boundary conditions
T_south =   0.
T_east  =  50.
T_north = 100.
T_west  =  75.0

T[:,1]   .= T_south
T[end,:] .= T_east
T[:,end] .= T_north
T[1,:]   .= T_west
# Note: I may want to fix the corners for plotting, but it's not needed for the solve.

# Loop until convergence
iter = 0
tol = 1e-6
iterMax = nx*ny*10

# since variables inside a loop are local, then we need to use this 'let' stuff to pass in the stuff from before.
let flag=0,tol=tol, iterMax=nx*ny*10, T=T
    while flag == 0
        # iter is 'global' since I want to reference it after the while loop ends.  This is needed if this while loop was inside a function.
        global iter += 1

        for i = 2:nx+1
            for j = 2:ny+1
                T[i,j] = 1/4*(T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1])
            end
        end

        # Compute R = Ax - b
        resid = zeros(nx+2,ny+2)
        for i = 2:nx+1
            for j = 2:ny+1
                resid[i,j] = -T[i,j] + 1/4*(T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1])
            end
        end

        # Since resid is a matrix, I'll use a norm
        if norm(resid) <= tol
            flag = 1
            println(iter)
        elseif iter >= iterMax
            flag = -1
        end

    end
end

@printf("iterations = %d\n",iter)
@printf("T = ")
display(T)
