function delta = data_adaptive_delta(x_mesh, ux_val, fx_val, threshold_u, threshold_f)
%% ux_val, fx_val: N * x_num matrices
% delta = support of hypothesis space H
sz = size(ux_val);
if numel(sz) == 2
    N = sz(1);
    T = 1;
    ux_val = reshape(ux_val, N, 1, sz(2));
    fx_val = reshape(fx_val, N, 1, sz(2));
else
    N = sz(1);
    T = sz(2);
end
supp_u = zeros(N, T, 2); supp_f = zeros(N, T, 2);
for n = 1:N
    parfor t = 1:T
        supp_u(n,t,:) = [x_mesh(find(abs(ux_val(n,t,:))>threshold_u, 1, 'first')), x_mesh(find(abs(ux_val(n,t,:))>threshold_u, 1, 'last'))];
        supp_f(n,t,:) = [x_mesh(find(abs(fx_val(n,t,:))>threshold_f, 1, 'first')), x_mesh(find(abs(fx_val(n,t,:))>threshold_f, 1, 'last'))];
    end
end
diff_supp = zeros(size(ux_val));
diff_supp(:,:,1) = supp_u(:,:,1) - supp_f(:,:,1); 
diff_supp(:,:,2) = supp_f(:,:,2) - supp_u(:,:,2);
[delta, ~] = max(diff_supp, [], 'all', 'linear');
delta = 1.1*delta;
end