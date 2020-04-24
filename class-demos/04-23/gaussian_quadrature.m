function [xi,w] = gaussian_quadrature(n_points)
% returns the natural points $\xi \in [-1,1]$ and the associated weights of
% 1D Gaussian quadrature

switch n_points
    
    case 1
        xi = 0;
        w  = 2;
        
    case 2
        xi = [-1/sqrt(3); 1/sqrt(3)];
        w  = [1; 1];
        
    case 3
        xi = [-3/sqrt(5); 0; 3/sqrt(5)];
        w = [5/9; 8/9; 5/9];        
        
        
    otherwise
        error('Order <%d> not implemented', d);
end