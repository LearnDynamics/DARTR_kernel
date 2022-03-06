function [basis_funs,knots] = Basis_bspline(degree,knots,knotsAug)
% compute the Bspline basis functions 
if knotsAug == 0.5 && knots(1)~=knots(2)
    knots = [knots(1)*ones(degree,1); knots(2:end)];
elseif knotsAug == 1 && knots(1)~=knots(2)&& knots(end)~=knots(end-1)
    knots = [knots(1)*ones(degree,1); knots(2:end-1);knots(end)*ones(degree,1)];
end

basis_funs  = localBasisFn_bspline(knots,degree); % B_i,n, i = 0, ... , k-degree-1, n = degree

%basisInfo.basis_funs = basis_funs; 
%basisInfo.degree     = degree;
%basisInfo.knots = knots; 
end


function basis_funs = localBasisFn_bspline(t,degree)

num_e = length(t);
basis_funs = cell(1, num_e - degree); % B_{i,degree}, i = 1, ... , k-degree-1

for i = 1:num_e - degree 
    basis_funs{i} = @(x) bspline_basis(i, degree, t, x);  
end
end

function y = bspline_basis(i, j, t, x)
% B-spline basis function value B(i,j) at x with konts t.
% i is the knots position: 1, 2, 3, ......
% j is the degree: 0, 1, 2, ......
y = zeros(size(x));

if j > 1
    b = bspline_basis(i,j-1,t,x);
    dn = x - t(i);
    dd = t(i+j -1) - t(i);
    if dd ~= 0  % indeterminate forms 0/0 are deemed to be zero
        y = y + b.*(dn./dd);
    end
    b = bspline_basis(i+1,j-1,t,x);
    dn = t(i+j) - x;
    dd = t(i+j) - t(i+1);
    if dd ~= 0
        y = y + b.*(dn./dd);
    end
else
    y(t(i) < x & x<= t(i + 1) ) = 1;
end
end

