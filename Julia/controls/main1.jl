using Pkg
Pkg.activate(".")
using ControlSystems
using LinearAlgebra
import Polynomials
import DSP.conv

G = tf( 30, conv([1,10], [1,0.5, 3]) )

Ts = 2
OS = 50/100

A = [0 1 0;
     0 0 1;
     -30 -8 -10.5]
B = [0.; 0; 1]
C = [30. 0 0]
D = 0.

ζ = sqrt( log(OS)^2/(pi^2 + log(OS)^2))
ωn = 4/ζ/Ts

sx1 = -ζ*ωn + 1im*ωn*sqrt(1-ζ^2)
sx2 = conj(sx1)
sx3 = 20*real(sx1)

d = Polynomials.poly([sx1,sx2,sx3])
k = d.a[1:end-1] - G.matrix[1].den.a[1:end-1] |> real


T = ss( A-B*k', B, C, D )

stepplot(T)
