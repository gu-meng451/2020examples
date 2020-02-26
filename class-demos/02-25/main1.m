%% First Ex

clear all
close all
clc

%% System Parameters

l1 = 0.5; %[m]
l2 = 0.5;%[m]
m1 = 1; %[kg]
m2 = 1; %[kg]
g  = 9.81; %[m/s^2]

tf = 20; %[s]

%% Initial conditions
theta1_0 = 45*pi/180; %[rad]
theta2_0 = 60*pi/180; %[rad]

theta1d_0 = 10; %[rad/s] 
theta2d_0 = -10; %[rad/s]

%%
% Let's transform these coordinate-based IC to our states
p1_0 = l1.^2.*m1.*theta1d_0 + m2.*(2*l1.^2.*theta1d_0 + 2*l1.*l2.*cos(theta1_0 - theta2_0).*theta2d_0)/2;
p2_0 = m2.*(2*l1.*l2.*cos(theta1_0 - theta2_0).*theta1d_0 + 2*l2.^2.*theta2d_0)/2;

z0 = [theta1_0; theta2_0; p1_0; p2_0 ];


%% Let's solve
sol = ode45( @(t,z) doublepen(t,z,l1,l2,m1,m2,g), [0,tf], z0);

%%
% plot(sol.x, sol.y(1,:))

%%
th1 = @(t) deval( sol, t, 1 );
th2 = @(t) deval(sol, t, 2);

%%
fig = figure();
ax = axes('Parent', fig);
hold(ax,'on');
set(ax, 'DataAspectRatio', [1,1,1]);
xlabel(ax, 'x')
ylabel(ax, 'y')
xlim(ax,[-1,1]*(l1+l2) );
ylim(ax,[-1,1]*(l1+l2) );

x1 = @(t) l1*sin(th1(t));
y1 = @(t) -l1*cos(th1(t));

x2 = @(t) x1(t) + l2*sin(th2(t));
y2 = @(t) y1(t) - l2*cos(th2(t));


h1 = plot(ax, [0, x1(0)], [0, y1(0)], 'LineWidth', 3, 'color', [200, 16, 46]/255 );
h2 = plot(ax, [x1(0), x2(0)], [y1(0), y2(0)], 'LineWidth', 3, 'color', [4,30,66]/255 );

h_title = title(ax, 'place holder');

nframes = 50*tf;
time = linspace(0, tf, nframes);
for i = 1:nframes
   
    h1.XData = [0, x1(time(i))];
    h1.YData = [0, y1(time(i))];
    
    h2.XData = [x1(time(i)), x2(time(i))];
    h2.YData = [y1(time(i)), y2(time(i))];
    
    h_title.String = sprintf('Time t = %.3f [s]', time(i) );
    
    pause(0.1)
end

