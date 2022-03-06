function [basis_funs,basisInfo] = Basis_fourier(K, odd_even, supp)
% compute the fourier basis functions

%xmid = 0; 
dx   = 2 * pi; 
omega = 2*pi/dx; 
basis_funs = cell(K,1);
dbasis_funs = cell(K,1);
ddbasis_funs = cell(K,1);
s       = floor(K/2);
if exist('supp', 'var') && (supp(end)-supp(1) > 0) 
    d = 0.01; % smooth range
    % smooth indicator fun
%     A1 = [pi^3 pi^2 pi 1; (pi-d)^3 (pi-d)^2 (pi-d) 1; 3*pi^2 2*pi 1 0;3*(pi-d)^2 2*(pi-d) 1 0];
%     cubic_coef1 = lsqlin(A1, [0; 1; 0; 0]);
%     a3 = cubic_coef1(1); a2 = cubic_coef1(2);a1 = cubic_coef1(3);a0 = cubic_coef1(4);
%     cubic_right = @(x) a3*x.^3 +a2*x.^2 +a1*x +a0; 
%     dcubic_right = @(x) 3*a3*x.^2 + 2*a2*x +a1;
%     ddcubic_right = @(x) 6*a3*x + 2*a2;
%     A2 = [(-pi)^3 (-pi)^2 (-pi) 1; (-pi+d)^3 (-pi+d)^2 (-pi+d) 1; 3*(-pi)^2 2*(-pi) 1 0;3*(-pi+d)^2 2*(-pi+d) 1 0];
%     cubic_coef2 = lsqlin(A2, [0; 1; 0; 0]);
%     a3 = cubic_coef2(1); a2 = cubic_coef2(2);a1 = cubic_coef2(3);a0 = cubic_coef2(4);
%     cubic_left = @(x) a3*x.^3 +a2*x.^2 +a1*x +a0; 
%     dcubic_left = @(x) 3*a3*x.^2 + 2*a2*x +a1;
%     ddcubic_left = @(x) 6*a3*x + 2*a2;
%     indicator = @(x) (x>=supp(1)+d & x<=supp(end)-d) +  cubic_left(x).*(x>=supp(1) & x<=supp(1)+d)+cubic_right(x).*(x<=supp(end) & x>=supp(end)-d);
%     dindicator = @(x) dcubic_left(x).*(x>=supp(1) & x<=supp(1)+d)+dcubic_right(x).*(x<=supp(end) & x>=supp(end)-d);
%     ddindicator = @(x) ddcubic_left(x).*(x>=supp(1) & x<=supp(1)+d)+ddcubic_right(x).*(x<=supp(end) & x>=supp(end)-d); 
    % piecewisesmooth indicator fun
%     indicator = @(x) (x>supp(1)+d & x<supp(end)-d) + (x-supp(1)).*(x>=supp(1) & x<=supp(1)+d)/d +(supp(end)-x).*(x<=supp(end) & x>=supp(end)-d)/d;
%     dindicator = @(x)  (x>=supp(1) & x<=supp(1)+d)/d - (x<=supp(end) & x>=supp(end)-d)/d;
%     ddindicator = @(x) zeros(size(x));
    % discontinuous indicator fun
    indicator = @(x) (x>supp(1) & x<supp(end));
    dindicator = @(x) zeros(size(x)) ;
    ddindicator = @(x) zeros(size(x));

else
    indicator = @(x) ones(size(x));
    dindicator = @(x) zeros(size(x));
end
switch odd_even
    case 'sincos'
        for i = 1:s
            basis_funs{i,1} = @(x) cos(i*x*omega).*indicator(x);
            dbasis_funs{i,1} = @(x) -i*omega*sin(i*x*omega).*indicator(x) + cos(i*x*omega).*dindicator(x);
            ddbasis_funs{i,1} = @(x) -i*omega*i*omega*cos(i*x*omega).*indicator(x) -i*omega*sin(i*x*omega).*dindicator(x) -i*omega*sin(i*x*omega).*dindicator(x)+ cos(i*x*omega).*ddindicator(x);

        end
        for i = s+1:K
            j = i-s;
            basis_funs{i,1} = @(x) sin(omega*j*x).*indicator(x);
            dbasis_funs{i,1} = @(x) j*omega*cos(j*x*omega).*indicator(x)+sin(omega*j*x).*dindicator(x);
            ddbasis_funs{i,1} = @(x) -j^2*omega^2*sin(j*x*omega).*indicator(x)+j*omega*cos(j*x*omega).*dindicator(x)+j*omega*cos(omega*j*x).*dindicator(x)+sin(omega*j*x).*ddindicator(x);
        end
    case 'sine'
        for j = 1:K
             basis_funs{j,1} = @(x) sin(omega*j*x).*indicator(x);
            dbasis_funs{j,1} = @(x) j*omega*cos(j*x*omega).*indicator(x)+sin(omega*j*x).*dindicator(x);
            ddbasis_funs{j,1} = @(x) -j^2*omega^2*sin(j*x*omega).*indicator(x)+j*omega*cos(j*x*omega).*dindicator(x)+j*omega*cos(omega*j*x).*dindicator(x)+sin(omega*j*x).*ddindicator(x);
        
        end
    case 'cosine'
        for j = 1:K
            i = j;%j - 1;
            basis_funs{i,1} = @(x) cos(i*x*omega).*indicator(x);
            dbasis_funs{i,1} = @(x) -i*omega*sin(i*x*omega).*indicator(x) + cos(i*x*omega).*dindicator(x);
            ddbasis_funs{i,1} = @(x) -i*omega*i*omega*cos(i*x*omega).*indicator(x) -i*omega*sin(i*x*omega).*dindicator(x) -i*omega*sin(i*x*omega).*dindicator(x)+ cos(i*x*omega).*ddindicator(x);
        end
end
  
basisInfo.basis_funs = basis_funs; 
basisInfo.type       = 'Fourier';  
basisInfo.odd_even   = odd_even;
basisInfo.dbasis_funs = dbasis_funs; 
basisInfo.ddbasis_funs = ddbasis_funs; 
end


