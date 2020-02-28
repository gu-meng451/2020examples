%% Constrained Lag
clear all
close all
clc

%%
g = 9.81; %[m/s^2]
l = 1; %[m]
tf = 10;

%% Initial conditions
theta0 = 45*pi/180;

x0 = l*sin(theta0);
y0 = -l*cos(theta0);
dx0 = 0;
dy0 = 0;


%%
% z= [x,y,xd, yd]
z0= [x0, y0, dx0, dy0]';
F = @(t,z) [z(3);
            z(4);
            [-z(2), z(1); z(1), z(2)]\[-g*z(1); -z(3)^2 - z(4)^2] ];

sol = ode45(F, [0,tf], z0);

%%
plot( sol.x, sol.y,'-', 'marker', '.')

%%
ll= @(t) deval(sol, t, 1).^2 + deval(sol, t, 2).^2 - l^2;
figure
fplot( ll, [0, tf])