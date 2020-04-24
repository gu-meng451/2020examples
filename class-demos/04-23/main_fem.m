%% finite element example

%%
clear all
close all
clc

%% define parameters
L = 1; % [m]
A = 0.002581; %[m^2] ~ (2in)^2
E = 69e9; % [Pa]

famp = 2000/E/A;
f = @(x) (1-x/L).^2*famp; % [N]

%% Boundary conditions
% $x = 0$, Dirchlet at the wall
u0 = 0;

%%
% $x = L$ Neumann at the free end
h = 0/E/A;


%% Make the mesh
n = 10;
x = linspace(0,L,n);

%% Global stiffness matrix
K = zeros( n, n);

% number of elements:
nel = n-1;

for e = 1:nel
    j = e:e+1;
    A = e;
    ke = 1/(x(A+1) - x(A))*[1, -1; -1, 1];
    K(j,j) = K(j,j) + ke;
end


%% build Right-hand-side (rhs)
R = zeros(n,1);
quad_n = 3;
[quad_xi,quad_w] = gaussian_quadrature(quad_n);
N = @(xi) [(1-xi)/2, (1+xi)/2];

for e = 1:nel
    j = e:e+1;
    A = e;
    xe = [x(A); x(A+1)];   
    J = (x(A+1) - x(A))/2;
    fe = zeros(2,1);

    for l = 1:quad_n
        xi = quad_xi(l);
        fe = fe + N(xi).'*f( N(xi)*xe )*J*quad_w(l);
    end
    
    R(j) = R(j) + fe;
    
end

%% Apply Dirchlet Boundary Conditions
% since u(A=1) is the only Dirchlet BC, we can simple move that
% contribution to the RHS
K_ff = K(2:n, 2:n);
K_rf = K(2:n, 1);
R_ff = R(2:n);

u_f = K_ff\(R_ff - K_rf*u0);

% rebuild u
u = [u0; u_f];

%% plot results
plot( x, u, 'o-', 'DisplayName', 'FEM Solution', 'LineWidth', 2 )
hold on

%% analytical solution
F = @(x) -famp*(x.^4/12/L^2 - x.^3/3/L + x.^2/2);
c2 = u0;
F1 = @(x) -famp*(x.^3/3/L^2 - x.^2/L + x);
c1 = h - F1(L);
sol = @(x) F(x) + c1*x + c2;
plot( x, sol(x), '--', 'LineWidth', 2, 'DisplayName', 'Analytic Solution')

%%
legend('show', 'location', 'SouthEast');