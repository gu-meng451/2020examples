function dz = myode(t,z,A,B,C,K,L,r)

n = size(A,1);

x = z(1:n);
y = C*x;
xhat = z(n+1:2*n);

%% Build Observer
u = r(t) - K*xhat;
yhat = C*xhat;
dxhat = A*xhat + B*u + L*(y - yhat);

%% The 'actual system' 
% We don't know x

dx = A*x + B*u;

%% Pack up the output
dz = [dx; dxhat];