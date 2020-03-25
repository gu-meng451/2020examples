%function x = ForwardEuler(f, tf, h, x0)


clear all
close all
clc

%%
f = @(t,x) -2*x;

x0 = 1;


%%
h = 0.1/2/2;
tf = 10;


%%
time = 0:h:tf;

%%
n = length(time);
x = nan( n,1 );

%%
x(1) = x0;

for i = 1:n-1
   x(i+1) = x(i) + h*f(time(i), x(i) ); 
end

%%
plot(time, x, 'DisplayName', 'Forward Euler')
hold on
fplot(@(t) exp(-2*t)*x0, [0,tf], '--', ...
    'DisplayName', 'Exact Solution')
legend('on')






% end