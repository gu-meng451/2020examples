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