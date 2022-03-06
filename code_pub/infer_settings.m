function [inferInfo] = infer_settings(obsInfo, deg, M,support_r)
% 
x_mesh = obsInfo.x_mesh; 
dx     = x_mesh(2) - x_mesh(1);
supp_H = support_r;
% supp_H =  obsInfo.supp_H;
inferInfo.Bmat_type = 'Lebesgue'; % or 'L2rho', or 'Unifrho'    % it may be overwritten later, since infer setting is not used for vector estimator
                          % we could save all Bmat's here. Now re-computed at learn_kernel2 and saved in Est
%% H = span{phi_1, ..., phi_M} on supprt [dx, xmax - u_supp(end) ]
%M = 10;
%supp_H = [0 max([default_x_mesh(ind_max) - u_supp(2), u_supp(1) - default_x_mesh(ind_min)])];
if ~exist('M', 'var')
    M = min(ceil(supp_H(2)/dx)*2, 300);
end
basis_Type = 'Bspline'; 
switch basis_Type
    case {'Bspline'}
        degree   = deg; % 1, 2, 3...
        knotsAug = 0.5; 
        %dist_x = abs(x_mesh' - x_mesh); dist_x = reshape(dist_x,1,[]);
        %dist_x(dist_x>supp_H(end)) = [];
        %quantile_grid = quantile(dist_x, M+degree-4);
        %knots = [supp_H(1) quantile_grid supp_H(end)]';
        knots = linspace(supp_H(1), supp_H(end), M+degree-2)';    %% --- not matching vector estimator  
        [basis_funs, knots] = Basis_bspline(degree,knots,knotsAug);
        inferInfo.degree = degree;
        inferInfo.basis_funs = basis_funs;
        inferInfo.basis_funs_knots = knots;
        inferInfo.basis_str   = [basis_Type,'_','quantileKnots',num2str(knots(1)), '_',num2str(floor(knots(end))), '_M_',num2str(numel(basis_funs)),'_deg_',num2str(degree)];
        inferInfo.basis_str   = [basis_Type,'_','UnifKnots',num2str(knots(1)), '_',num2str(floor(knots(end))), '_M_',num2str(numel(basis_funs)),'_deg_',num2str(degree)];

end
inferInfo.supp_H = supp_H;
inferInfo.delta  = obsInfo.delta;
inferInfo.B_Lebesgue = 0; % computeB_Lebesgue(basis_funs,supp_H,knots, degree); 

end


function B_Lebesgue = computeB_Lebesgue(basis_funs,supp_H, basis_knots, basis_deg)
% compute Bmat with Lebesgue measure 
% B_Lebesgue(i, j) = <phi_i, phi_j>_L2
nBasis = numel(basis_funs);
Bspline_flag = exist('basis_knots', 'var') && exist('basis_deg', 'var');
B_Lebesgue = zeros(nBasis);
parfor i = 1:nBasis
    for j = 1:nBasis
        if Bspline_flag
            xmin = min([basis_knots(i), basis_knots(j)]);
            xmax = max([basis_knots(i+basis_deg), basis_knots(j+basis_deg)]);
        else
            xmin = supp_H(1); xmax = supp_H(2);
        end
        fun = @(x) basis_funs{i}(x) .* basis_funs{j}(x);
        B_Lebesgue(i,j) = integral(fun, xmin, xmax);   % upper triangular       
    end
end
%B_Lebesgue = triu(B_Lebesgue) + tril(B_Lebesgue', -1); % symmetric
end