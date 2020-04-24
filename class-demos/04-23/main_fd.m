%% Example: finite difference bar

%%
clear all
close all
clc

%% define parameters
L = 1; % [m]
A = 0.002581; %[m^2] ~ (2in)^2
E = 69e9; % [Pa]

famp = 2000;
f = @(x) (1-x/L).^2*famp; % [N]

%% Boundary conditions
% $x = 0$, Dirchlet at the wall
u0 = 0;

%%
% $x = L$ Neumann at the free end
h = 0/E/A;

%%

%% Define the grid
n_list = [ 10, 20, 40];

for k = 1:length(n_list)

    n = n_list(k);
    x = linspace(0,L,n);
    dx = x(2) - x(1);
    
    %% Solve
    % I'll use Gauss-Siedel because it's easy
    % initialize u
    if k == 1
        u = zeros(n,1);
    else
        % use old solutiont to speed things up
        u = interp1( linspace(0,L, n_list(k-1) ), u, x);
    end
    
    % apply Dirchlet BC
    u(1) = u0;
    
    % iteratively solve
    tol = 1e-12;
    iter_max = 3*n^3;
    iter = 0;
    flag = 0;
    
    res = zeros(n,1);
    
    while flag == 0
        
        iter = iter+1;
        
        for i = 2:n-1
            u(i) = ( u(i+1) + u(i-1) + dx^2/E/A*f(x(i)) )/2;
        end
        u(n) = u(n-1);
        
        % check residuals
        for i = 2:n-1
            res(i) = ( u(i+1) + u(i-1) + dx^2/E/A*f(x(i)) )/2 - u(i);
        end
        res(n) = u(n-1) - u(n);
        
        if norm(res)/norm(u) < tol
            flag = 1;
            
        elseif iter >= iter_max
            flag = -1;
        end
        
    end
    
    %%
    plot( x, u, 'o-', 'DisplayName', sprintf('n = %4d',n), ...
        'LineWidth', 2)
    hold on
    drawnow
end

%% compute the analytic solution
F = @(x) -famp*(x.^4/12/L^2 - x.^3/3/L + x.^2/2)/E/A;
c2 = u0;
F1 = @(x) -famp*(x.^3/3/L^2 - x.^2/L + x)/E/A;
c1 = h - F1(L);
sol = @(x) F(x) + c1*x + c2;
plot( x, sol(x), '--', 'LineWidth', 2, 'DisplayName', 'Analytic Solution')

%%
legend('show', 'location', 'SouthEast');