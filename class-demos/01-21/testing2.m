%%
clear all
close all
clc

%%
x = -70.3
y_true = exp(x);

[y1,i1] = Exp1(x);
i1

err1 = ( y1 - y_true )/y_true

%%
[y2,i2] = Exp2(x);
i2

err2 = ( y2 - y_true)/y_true

%% Timing
X = 10*(2*rand(10000,1)-1);

tic
for i = 1:length(X)
    y = Exp2(X(i));
end
t1 = toc;
% Let's convert from sec to nsec
fprintf( 'Exp2 avg time: %6.1f [ns]\n', (t1*1e9)/length(X) );

tic
for i = 1:length(X)
    y = exp(X(i));
end
t2 = toc;
fprintf( 'exp avg time:  %6.1f [ns]\n', (t2*1e9)/length(X) );


