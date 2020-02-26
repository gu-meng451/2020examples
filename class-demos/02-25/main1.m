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


%% Initial conditions
theta1_0 = 45*pi/180; %[rad]
theta2_0 = 60*pi/180; %[rad]

theta1d_0 = 10; %[rad/s] 
theta2d_0 = -10; %[rad/s]

%%
% Let's transform these coordinate-based IC to our states
p1_0 = l1.^2.*m1.*theta1d_0 + m2.*(2*l1.^2.*theta1d_0 + 2*l1.*l2.*cos(theta1_0 - theta2_0).*theta2d_0)/2;
p2_0 = m2.*(2*l1.*l2.*cos(theta1_0 - theta2_0).*theta1d_0 + 2*l2.^2.*theta2d_0)/2;


