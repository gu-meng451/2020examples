%% Constrained Lagrange's Equation
clear all
close all
clc

%% Define paramters
g = 9.81; %[m/s^2]
l = 1; %[m]
tf = 10; %[s]

%% Initial conditions
theta0 = 45*pi/180;

%%
% in terms of the generalized coordinates (x,y)
x0 = l*sin(theta0);
y0 = -l*cos(theta0);
dx0 = 0;
dy0 = 0;

%% Solving the ODE IVP:
% z= [x,y,xd, yd]
z0= [x0, y0, dx0, dy0]';
F = @(t,z) [z(3);
            z(4);
            [-z(2), z(1); z(1), z(2)]\[-g*z(1); -z(3)^2 - z(4)^2] ];

sol = ode45(F, [0,tf], z0);

%%
% Simple plot to see the solutions
plot( sol.x, sol.y,'-', 'marker', '.')
xlabel('Time t [s]')
ylabel('States')
legend(["x", "y", "x dot", "y dot"])

%% What about the constraint?
h = @(t) (deval(sol, t, 1).^2 + deval(sol, t, 2).^2 - l^2)/l^2;
fig = figure();
ax = axes('Parent', fig);
fplot(ax, h, [0, tf], 'DisplayName', 'ode45');
xlabel(ax, 'Time t [s]');
ylabel(ax, 'Normalized Constraint h/l^2 = (x^2 + y^2 - l^2)/l^2');
hold(ax, 'on');

%%
% What if we used other solvers?
vsol = @(ode) ode( F,[0,tf], z0);
solver_list  = {@ode23t, @ode23tb, @ode113};
solver_names = ["ode23t", "ode23tb", "ode113"];

for i = 1:length(solver_list)
   sol = vsol(solver_list{i});
   h = @(t) (deval(sol, t, 1).^2 + deval(sol, t, 2).^2 - l^2)/l^2;
   fplot(ax, h, [0,tf], 'DisplayName', solver_names(i) );
end
legend(ax, 'show')

%% What is the Hamiltonian doing?
H = @(z) 1/2*(z(3,:).^2 + z(4,:).^2) + g*z(2,:);

solver_list  = {@ode45, @ode23t, @ode23tb, @ode113};
solver_names = ["ode45", "ode23t", "ode23tb", "ode113"];

fig3= figure();
ax3 = axes('Parent', fig3);
hold(ax3, 'on');
for i = 1:length(solver_list)
   sol = vsol(solver_list{i});
   plot(ax3, sol.x, H(sol.y)/H(z0), ...
       'DisplayName', solver_names(i) );
end
legend(ax3, 'show')
xlabel(ax3, 'Time t [s]');
ylabel(ax3, 'Normalized Hamiltonian H(t)/H(0)');

