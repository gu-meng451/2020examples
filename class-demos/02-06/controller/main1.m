%%
clear all
close all
clc


%%
G = tf( 30, conv([1,10], [1,0.5, 3]) )
Ts = 2;
OS = 50/100;

%%
A = [0, 1, 0;
     0, 0, 1;
     -30, -8,-10.5];
B = [0; 0; 1];
C = [30, 0, 0];
D = 0;

%%

zeta = sqrt( log(OS)^2/(pi^2 + log(OS)^2));
wn = 4/zeta/Ts;

sx1 = -zeta*wn + 1j*wn*sqrt(1-zeta^2)
sx2 = conj(sx1)
sx3 = 20*real(sx1)

%%
d = flip( poly([sx1,sx2,sx3]) )
a = [30, 8,10.5]
k = d(1:end-1) - a

%%
T = ss( A-B*k, B, C, D )
step(T)
stepinfo(T)


%%
C = eye(3)
T = ss( A-B*k, B, C, D )
step(T)
stepinfo(T)

[x,t] = step(T);

figure
plot( t, 1-x*k')




