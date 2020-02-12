%% 11-Feb: Observer example
% A continuation from last class' example.
clear all
close all
clc

%% Starting parameters
% Transfer function and target performance parameters
G = tf( 30, conv([1,10], [1,0.5, 3]) )
Ts = 2;
OS = 50/100;

%% TF to SS
A = [0, 1, 0;
     0, 0, 1;
     -30, -8,-10.5];
B = [0; 0; 1];
C = [30, 0, 0];
D = 0;

%% CL Poles

zeta = sqrt( log(OS)^2/(pi^2 + log(OS)^2));
wn = 4/zeta/Ts;

sx1 = -zeta*wn + 1j*wn*sqrt(1-zeta^2)
sx2 = conj(sx1)
sx3 = 5*real(sx1)

%% Design the controller (already in controller canonical form)
d = flip( poly([sx1,sx2,sx3]) )
a = [30, 8,10.5]
k = d(1:end-1) - a

%% Check results
T = ss( A-B*k, B, C, D )
step(T)
stepinfo(T)


%% Observer Design
% Make poles of A-LC faster than A-BK:
so1 = 10*sx1;
so2 = 10*sx2;
so3 = 10*sx3;
do = flip( poly([so1,so2,so3]) )

% Cast the observer problem as the transpose of the controller problem
% so A'->A and C'->B
At = A';
Bt = C';

% build controllability matrix in terms of state's we're working in
CMx = [Bt, At*Bt, At^2*Bt]

coeffCharPoly = poly(eig(At)) % same as they were in G, but in a different problem we may not know that yet
a = flip(coeffCharPoly(2:end))

Abar = [0,1,0;
        0, 0, 1;
        -a]; % put the coeffs in the last row
Bbar = [0;0;1];

% build controllability matrix in terms of controller canonical form states
% z
CMz = [Bbar, Abar*Bbar, Abar^2*Bbar]
P = CMx/CMz

% design controller gains
Lzt = do(1:end-1) - a

% transform from states z and transpose to get correct dimension (L is a
% column)
Lx = (Lzt/P)'

%% Check
% using the built-in function `place` on the transpose of the controller
% problem
Lx_check = place(A',C',[so1,so2,so3])'


%% Let's test the observer:

% actual initial conditions in the system
x0 = [1; -1; 10];

% our starting guess for our observer
xhat0 = [0;0;0];

% what we want the system to do
r = @(t) 0;
% if we set r(t) = 1, why doesn't our output go to 1?

% pack up the total initial states
z0 = [ x0; xhat0];

% period of integration
tspan = [0, 2];

% solve the system using ode45, because we haven't learned why not to use
% it yet,
sol = ode45( @(t,z) myode(t,z,A,B,C,k,Lx,r), tspan, z0);


%% Plotting the results 
% Really rough first plot
figure
plot(sol.x, sol.y(1,:), 'LineWidth', 3 )
hold on
plot(sol.x, sol.y(3+1,:), '--', 'LineWidth', 3 )
xlabel('Time t')
ylabel('state 1')


%% Fancy plot
% Because plots matter
figure
n = size(A,1);

ax1 = subplot(4,1,1);
y = @(t) C*deval(sol, t, 1:n );
yhat = @(t) C*deval(sol, t, n+1:2*n );
fplot(ax1, y, [0, 2], 'LineWidth', 3, 'DisplayName','y')
hold(ax1, 'on')
fplot(ax1, yhat, [0, 2], '--','LineWidth', 3, 'DisplayName','yhat')
legend(ax1, 'show')

for i = 1:n
    ax(i) = subplot(4,1,i+1);
    x = @(t) deval(sol, t, i);
    xhat = @(t) deval(sol, t, n+i );
    fplot(ax(i), x, [0, 2], 'LineWidth', 3)
    hold(ax(i), 'on')
    fplot(ax(i), xhat, [0, 2], '--', 'LineWidth', 3)
    ylabel(sprintf('State x_%d', i))
end
xlabel(ax(end), 'Time t')
