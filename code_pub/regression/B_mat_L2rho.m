function BL2rho = B_mat_L2rho(inferInfo, rho, SAVE_DIR)
% Output: B(i, j) = <phi_i, phi_j>_{L^2}
basis       = inferInfo.basis_funs; 
basis_knots = inferInfo.basis_funs_knots; 
basis_deg   = inferInfo.degree;
basis_str   = inferInfo.basis_str;
nBasis      = numel(basis);
delta       = inferInfo.supp_H(end);
filename    = [SAVE_DIR, 'B_mat_', basis_str, '.mat'];
rho_val     = rho.rho_val; rSeq = rho.r_seq;   

% flag = 0; 
% if isfile(filename);    load(filename, 'BL2rho');
%    if ~exist('BL2rho','var'); flag =1; end 
% end

% if ~isfile(filename) || flag ==1 
BL2rho = zeros(nBasis);
for i = 1:nBasis
    for j = i:nBasis
        if contains(inferInfo.basis_str, 'Bspline')
            xmin = min([basis_knots(i), basis_knots(j)]);
            xmax = max([basis_knots(i+basis_deg), basis_knots(j+basis_deg)]);
        elseif contains(inferInfo.basis_str, 'Bernstein')
            xmin = 0;
            xmax = delta;
        end
        fun = @(x) basis{i}(x) .* basis{j}(x);
        %B(i,j) = inner_prod(basis{i}, basis{j}, [xmin xmax]);
        %  B(i,j) = integral(fun, xmin, xmax);
        
        fun_val= fun(rSeq);
        BL2rho(i,j) = fun_val*rho_val;
        
        BL2rho(j,i) = BL2rho(i,j);
    end
end
dr = rSeq(2)-rSeq(1);
BL2rho = BL2rho*dr;
if min(eig(BL2rho)) <1e-20; BL2rho = BL2rho+eye(nBasis)*1e-20; end
%     save(filename, 'BL2rho','-append');
% end


end


