%% Example: using Euler

clear all
close all
clc

%% Problem to solve
f = @(t,x) -2*x;
x0 = 1;

% known exact solution:
X = @(t) exp(-2*t)*x0;


%% Solve using Forward Euler:

%%
% Setup the time stuff
h = 0.1;
tf = 4;

time = 0:h:tf;

%% 
% Build the storage container (this time I want to save everything)
n = length(time);
x = nan( n,1 );

%%
% Initial step
x(1) = x0;

%%
% Step through to the end
for i = 1:n-1
   x(i+1) = x(i) + h*f(time(i), x(i) ); 
end

%% Plot the results

plot(time, x, 'DisplayName', 'Forward Euler')
hold on
fplot(@(t) exp(-2*t)*x0, [0,tf], '--', ...
    'DisplayName', 'Exact Solution')
legend('show')
xlabel('Time t [s]')
ylabel('x(t) [units]')

