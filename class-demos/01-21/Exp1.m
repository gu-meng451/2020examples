
function [y,i] = Exp1(x)
%%

% y = exp(x)


%%
y = 1;
newterm = 1;

i = 0;
flag = 0;
tol = 1e-12;
max_iter = 300;

while flag == 0
   i = i+1;
   
   newterm = newterm*x/i;
   
   y = y + newterm;
   
   rel_err = abs( newterm );
   
   if rel_err <= tol
       flag = 1;
       
   elseif i >= max_iter
       flag= -1;
       error('Failed to converge')
   end
    
end
