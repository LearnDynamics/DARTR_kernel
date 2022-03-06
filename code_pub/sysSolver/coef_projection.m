function [coef,xi,phiVal,trueFn_val,projFn_val] = coef_projection(bFn,trueFn, xgrid)
% find the best coefficients given a set of basis fucntions in L2(\mu)
%        trueFn(xi) = sum_k c_k phi_k(xi),    regression in L2 to get c_k
%        < trueFn, phi_l>_{L2} = <sum_k c_k phi_k, phi_l>_{L2}
xi = xgrid;
%xi = (xi(1:end-1) + xi(2:end) )/2;  % use mid-point for the L2 integral
xi = xi(1:end);  % use right-point for the L2 integral

xi = xi';
n  = length(bFn); 

A    = zeros(n,n); b = zeros(n,1);

trueFn_val = trueFn(xi);
phiVal     = zeros(length(xi),n);
for i=1:n
        phi_i   = bFn{i}; 
        phi_val = phi_i(xi);   phiVal(:,i) = phi_val; 
        b(i)    = sum(phi_val.* trueFn_val);  
    for j=i:n
        phi_j    = bFn{j}; 
        phi_j_val= phi_j(xi);
        A(i,j)   = sum(phi_val.*phi_j_val); 
        A(j,i)   = A(i,j);
    end
end
coef       = pinv(A)*b; 
projFn_val = phiVal*coef;


plotON = 0;
if plotON==1
    figure;
    plot(xi,projFn_val,'b-','linewidth',1);hold on;
    plot(xi,trueFn_val,'r:','linewidth',1);hold on;
    legend('Projection','True');
end
end