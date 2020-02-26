
function doublepen!(dz, z, p, t)
    # unpack the parameters
    l₁, l₂, m₁, m₂, g = p

    # unpack the states
    z1,z2,z3,z4 = z

    # compute the derivatives
    dz[:] = [(-l₁*z4*cos(z1 - z2) + l₂*z3)/(l₁^2*l₂*(m₁ - m₂*cos(z1 - z2)^2 + m₂)),
    (l₁*m₁*z4 + l₁*m₂*z4 - l₂*m₂*z3*cos(z1 - z2))/(l₁*l₂^2*m₂*(m₁ - m₂*cos(z1 - z2)^2 + m₂)),
     (-g*l₁^3*l₂^2*m₁^3*sin(z1) + 2*g*l₁^3*l₂^2*m₁^2*m₂*sin(z1)*cos(z1 - z2)^2 - 3*g*l₁^3*l₂^2*m₁^2*m₂*sin(z1) - g*l₁^3*l₂^2*m₁*m₂^2*sin(z1)*cos(z1 - z2)^4 + 4*g*l₁^3*l₂^2*m₁*m₂^2*sin(z1)*cos(z1 - z2)^2 - 3*g*l₁^3*l₂^2*m₁*m₂^2*sin(z1) - g*l₁^3*l₂^2*m₂^3*sin(z1)*cos(z1 - z2)^4 + 2*g*l₁^3*l₂^2*m₂^3*sin(z1)*cos(z1 - z2)^2 - g*l₁^3*l₂^2*m₂^3*sin(z1) + l₁^2*m₁*z4^2*sin(2*z1 - 2*z2)/2 + l₁^2*m₂*z4^2*sin(2*z1 - 2*z2)/2 - l₁*l₂*m₁*z3*z4*sin(z1 - z2) - l₁*l₂*m₂*z3*z4*sin(z1 - z2)*cos(z1 - z2)^2 - l₁*l₂*m₂*z3*z4*sin(z1 - z2) + l₂^2*m₂*z3^2*sin(2*z1 - 2*z2)/2)/(l₁^2*l₂^2*(m₁^2 - 2*m₁*m₂*cos(z1 - z2)^2 + 2*m₁*m₂ + m₂^2*cos(z1 - z2)^4 - 2*m₂^2*cos(z1 - z2)^2 + m₂^2)),
     (-g*l₁^2*l₂^3*m₁^2*m₂*sin(z2) + 2*g*l₁^2*l₂^3*m₁*m₂^2*sin(z2)*cos(z1 - z2)^2 - 2*g*l₁^2*l₂^3*m₁*m₂^2*sin(z2) - g*l₁^2*l₂^3*m₂^3*sin(z2)*cos(z1 - z2)^4 + 2*g*l₁^2*l₂^3*m₂^3*sin(z2)*cos(z1 - z2)^2 - g*l₁^2*l₂^3*m₂^3*sin(z2) - l₁^2*m₁*z4^2*sin(2*z1 - 2*z2)/2 - l₁^2*m₂*z4^2*sin(2*z1 - 2*z2)/2 + l₁*l₂*m₁*z3*z4*sin(z1 - z2) + l₁*l₂*m₂*z3*z4*sin(z1 - z2)*cos(z1 - z2)^2 + l₁*l₂*m₂*z3*z4*sin(z1 - z2) - l₂^2*m₂*z3^2*sin(2*z1 - 2*z2)/2)/(l₁^2*l₂^2*(m₁^2 - 2*m₁*m₂*cos(z1 - z2)^2 + 2*m₁*m₂ + m₂^2*cos(z1 - z2)^4 - 2*m₂^2*cos(z1 - z2)^2 + m₂^2))]
end

function hamiltonian(sol, params, t)
    # unpack the parameters
    l1, l2, m1, m2, g = params

    # unpack the states
    z1,z2,z3,z4 = sol(t)

    # Build the hamiltonian
    H = (-2*g*l1^2*l2^2*m2*(m1 - m2*cos(z1 - z2)^2 + m2)^2*(l1*m1*cos(z1) + l1*m2*cos(z1) + l2*m2*cos(z2)) + m1*m2*(l1*z4*cos(z1 - z2) - l2*z3)^2 + m2^2*(l1*z4*cos(z1 - z2) - l2*z3)^2 - 2*m2*(l1*z4*(m1 + m2) - l2*m2*z3*cos(z1 - z2))*(l1*z4*cos(z1 - z2) - l2*z3)*cos(z1 - z2) + (l1*z4*(m1 + m2) - l2*m2*z3*cos(z1 - z2))^2)/(2*l1^2*l2^2*m2*(m1 - m2*cos(z1 - z2)^2 + m2)^2)

    return H
end
