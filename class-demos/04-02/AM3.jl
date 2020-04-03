module AM3

using LinearAlgebra

export am3


function am3(f, tf, h, x0; tol = 1e-4, iterMax = 200 )
        time = 0:h:tf
        n = length(time)
        p = length(x0)

        x = zeros(n,p)

        # Initial Condition
        x[1,:] .= x0

        # 1 Euler Step (this is cheap and possibly bad)
        fn = f(x[1,:], time[1])
        y = x[1,:] + h*fn

        # B-Euler for 1 step
        flag = 0
        iter = 0
        while flag == 0
            iter += 1
            y = x[1,:] + h*f(y, time[2])

            residual = norm( y - x[1,:] - h*f(y, time[2]) )

            if residual <= tol
                flag = 1
            elseif iter >= iterMax
                flag = -1
                error("Error: failed to converge")
            end
        end

        x[2,:] = y
        # AB2/AM3 the rest
        for i = 2:n-1

            # AB2 Predictor
            y = x[i,:] + h/2*( 3*f(x[i,:], time[i]) - f(x[i-1,:], time[i-1]) )


            flag = 0
            iter = 0
            while flag == 0
                iter += 1
                ynew = x[i,:] + h*( 5/12*f(y, time[i+1]) + 2/3*f(x[i,:], time[i]) - 1/12*f(x[i-1,:], time[i-1])  )
                err= norm(y - ynew)
                y = ynew
                if err <= tol
                    flag = 1
                elseif iter >= iterMax
                    flag = -1
                    error("Failed to converge") #TODO add more fail info here
                end
            end
            x[i+1,:] = y

        end

        return x, time

end


end
