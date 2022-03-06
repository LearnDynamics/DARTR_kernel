function [EstAdaptive_inferInfo,  idx_nonzero_phi_conv]= estimation_adaptive_infer_settings(Est1, obsInfo, inferInfo)

x_mesh = obsInfo.x_mesh;        dx = obsInfo.x_mesh_dx;
phi    = inferInfo.basis_funs;   u = obsInfo.u; 
residual = Est1.residual;
[~, ind_max_resd] = max(abs(residual));
x_max = x_mesh(ind_max_resd); % the point that has maximal (f - f_hat)(x) 

% f_hat(x) = sum c_hat_i * phi_conv ===> want to "improve" nonzero/large c_hat_i * phi_conv
% phi_conv(i) = integral phi_i(|y-x|)( u(y) -u(x) ) dy
phi_conv = zeros(size(phi));
for i = 1:numel(phi)
    phi_conv(i) = sum(phi{i}(abs(x_mesh-x_max)).*(u{1}(x_mesh) -u{1}(x_max)))*dx;
end
phi_conv = Est1.c_hat'.* phi_conv;
H_degree = inferInfo.degree;
idx_nonzero_phi_conv = find(abs(phi_conv)>dx);
if numel(idx_nonzero_phi_conv) == 0
    EstAdaptive_inferInfo = inferInfo;
    return
end
%idx1 = idx_nonzero_phi_conv(1); idx2 = idx_nonzero_phi_conv(end);
knots = inferInfo.basis_funs_knots;
mid = 0.5*(nonzeros(knots(idx_nonzero_phi_conv) + knots(idx_nonzero_phi_conv+1)));
if numel(mid) == 0
    n = nonzeros(knots);
    mid = 0.5 * n(1);
end 
% half those intervals
knots = sort([knots; mid]);
[basis_funs_EstAdaptive, knots] = Basis_bspline(H_degree,knots,0);
EstAdaptive_inferInfo = inferInfo;
EstAdaptive_inferInfo.degree = H_degree;
EstAdaptive_inferInfo.basis_funs = basis_funs_EstAdaptive;
EstAdaptive_inferInfo.basis_funs_knots = knots;
EstAdaptive_inferInfo.basis_str = ['Bspline_EstAdaptiveKnots',num2str(knots(1)), '_',num2str(floor(knots(end))), '_M_',num2str(numel(basis_funs_EstAdaptive)),'_deg_',num2str(H_degree)];
EstAdaptive_inferInfo.supp_H = [knots(1) knots(end)];

end
