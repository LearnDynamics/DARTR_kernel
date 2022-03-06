
function regressionData = getData4regression(ux_val,div_u,fx_val,dx,obsInfo,bdry_width,Index_xi_inUse,r_seq,data_str, normalizeOn)
 % pre-processData for regression: 
 %             min |L_phi[u]- f|^2 = c'*Abar*c - 2c'*b + bnorm_sq  
%{ 
L_phi[u](x) = \int phi(|y|)*g[u](x,y) dy 
            = \sum_r phi(r)* [ g[u](x,x+r) + g[u](x,x-r) ] dr 
 1. g[u](x,y) =  u(x+y)          integral kernel: 1D, Laplace Transform 
 2. g[u](x,y) =  [u(x+y)-u(x)]   nonlocal kernel
 3. g[u](x,y) = div(u(x+y))u(x)  particle interaction --- 1D, div(u(x+y)) = u'(x+y) 
====== In the inverse-integral-operator
  <L_phi[u], L_psi[u]> = \int_x L_phi[u](x) L_psi[u](x) dx
                       = \sum_{r,s} phi(r)psi(s) \int_x convl_gu(r,x)convl_gu(s,x)dx drds
         convl_gu(r,x) = [ g[u](x,x+r) + g[u](x,x-r) ]          ***** to be computed here
                G(r,s) = \int_x convl_gu(r,x)convl_gu(s,x)dx    *****
 - for function learning: we need the
                A(i,j) = <L_{phi_i}[u], L_{phi_j}[u]> = sum_{r,s} phi_i(r)phi_j(s) G(r,s) drds
                b(i)   = <L_{phi_i}[u], f> = int_{r} phi_i(r) sum_x convl_gu(r,x)f(x)dx  dr
                       >>>>>                                  *** convl_gu_f(r)   *****
 
 - for vector learning: we only need [with phi_i being the indicator at r_i ]
                A(i,j) = G(r_i,r_j) *dr*dr
                b(i)   = G*phi  = convl_gu_f(r_i) *dr
===>>> ***** in either case, we only need those *****
- compute  convl_gu(r,x) = [ g[u](x,x+r) + g[u](x,x-r) ] -- used in G, rho
- compute  G(r,s)    --- used in vector & function learning

- compute  rho(r)    --- used in all:  exploration measure rho(r) = \int   |convl_gu(r,x)| dx
%}


%% 1 load and normalize process data (u,f)
if numel(size(ux_val)) == 3               % u(n,t,x) 
    [case_num, T, x_num] = size(ux_val);   N = case_num*T;
    ux_val      = reshape(ux_val, case_num*T, x_num);
    fx_val  = reshape(fx_val, case_num*T, x_num);
elseif numel(size(ux_val)) == 2           % u(n,x) 
    [case_num, x_num] = size(ux_val);     N = case_num; T = 1;
end
% for nonLinearOpt: 
du = zeros(size(ux_val)); ddu = zeros(size(ux_val));
if strcmp(obsInfo.example_type, 'nonlinearOpt') || strcmp(obsInfo.example_type, 'mfOpt')
    x_mesh = obsInfo.x_mesh;
    u_supp = obsInfo.u_supp;
    dx = obsInfo.x_mesh_dx;
    %1: approx div_u, lap_u by central diff
    du(:, 2:end-1) = (ux_val(:,3:1:end)-ux_val(:,1:1:end-2))/2/dx; 
    ddu(:, 2:end-1) = (ux_val(:,3:1:end) - 2*ux_val(:,2:1:end-1) + ux_val(:,1:1:end-2))/dx/dx; 
    drop_idx = (x_mesh>(u_supp(end)-dx)&x_mesh<u_supp(end)+dx)|(x_mesh<(u_supp(1)+dx)&x_mesh>u_supp(1)-dx);
    ddu(:,drop_idx) = 0;du(:, drop_idx) = 0;% we are droping the data near boundary of u_supp
    ux_val(:, drop_idx) = 0;fx_val(:, drop_idx) = 0;
    %2: exact div_u, lap_u
%     du = div_u;
%     ddu = lap_u;
%     drop_idx = (x_mesh>(u_supp(end)-dx)&x_mesh<u_supp(end))|(x_mesh<(u_supp(1)+dx)&x_mesh>u_supp(1)); 
%     ddu(:,drop_idx) = 0;du(:, drop_idx) = 0;
%     ux_val(:, drop_idx) = 0;fx_val(:, drop_idx) = 0;

    % substract ddu from fx_val
    fx_val  = fx_val-ddu;
    %normalizeOn = 0;
end


% normalize by L2 norm if normalizeOn 
if normalizeOn 
    [ux_val, fx_val] = normalize_byL2norm(ux_val, fx_val, dx); 
    data_str   = [data_str,'NormalizeL2']; 
end
N_xi_inUse  = x_num-2*bdry_width;
N_r_seq     = length(r_seq);     % r_seq       = dx*(1:bdry_width);
if (length(Index_xi_inUse) ~= N_xi_inUse || (N_r_seq ~= bdry_width && bdry_width >0)) && ~contains(data_str, 'Particles')
    error('Index_xi_inUse does not match with data and bdry_width.  \n'); % error and terminate 
end
%% 2. get data for  regression:   
if bdry_width == 0
% 2.0. get data for classical regression:   
    bnorm_sq = sum(fx_val.^2);
    regressionData.rho_val = ones(N_r_seq,1)/N_r_seq/dx;   % exploration measure: uniform 
    regressionData.Gbar    = eye(N_r_seq);      
    regressionData.gu_f = fx_val';      
    regressionData.gu_f2= fx_val(1:2:end)';    % a downsampled estimator of gu_fN 
	regressionData.bnorm_sq   = bnorm_sq/N;
	regressionData.bdry_width = bdry_width; 
	regressionData.r_seq      = r_seq;  
	regressionData.data_str   = data_str;
elseif contains(data_str, 'Particles')
%2.1. partical system
	[N, pN, d] = size(ux_val);
    pair_dif = zeros(N, pN, pN, d);
    for i = 1:pN
        pair_dif(:, i, :, :) = ux_val(:, i, :) - ux_val;
    end
    pair_dist = sqrt(sum(pair_dif.^2, 4));
    
    pair_dist_flat = reshape(pair_dist, 1, []);
    pair_dist_flat(pair_dist_flat == 0) = [];
    
    L = max(pair_dist_flat);
    x_mesh = 0:dx:L;
    regressionData.pair_dist_flat = pair_dist_flat;
    regressionData.pair_dist = pair_dist;
    regressionData.pair_dif = pair_dif;
    regressionData.ux_val = ux_val;
    regressionData.fx_val = fx_val;
    regressionData.particles = obsInfo.particles;

else
% 2.2. get data for regression: convl_gu(r,x) = g[u](x,x+r) + g[u](x,x-r); 
    fun_g_vec = obsInfo.fun_g_vec;
    ind_p = 1:bdry_width;   ind_m = -(1:bdry_width);       % Index plus r  

    rhoN     = zeros(N_r_seq,N);                % N copies: exploration measure rho: rho(r) = \int   |convl_gu(r,x)| dx
    GbarN    = zeros(N_r_seq,N_r_seq,N);        % N copies: G(r,s) = \int_x convl_gu(r,x)convl_gu(s,x)dx 
    gu_fN    = zeros(N_r_seq,N);                % N copies:  sum_x convl_gu(r,x)f(x)dx
    gu_fN2   = gu_fN;                           % a downsampled estimator of gu_fN 
    bnorm_sq = 0; 
    for nn=1:N                                  % compute these terms for each u(x) 
         u1  = ux_val(nn,:); f1 = fx_val(nn,Index_xi_inUse)';
         du1 = du(nn,:);
         convl_gu = zeros(N_r_seq,N_xi_inUse);  val_abs = convl_gu;  
         for k  = 1:N_xi_inUse
             temp_p        = fun_g_vec(u1,du1, k+bdry_width, ind_p); 
             temp_m        = fun_g_vec(u1,du1, k+bdry_width, ind_m); 
             convl_gu(:,k)   = (temp_p + temp_m)';
             val_abs(:,k)  = (abs(temp_p) + abs(temp_m))';     
         end
         rhoN(:,nn)     = mean(val_abs,2); 
         GbarN(:,:,nn)  = convl_gu* convl_gu'; 
         gu_fN(:,nn)    = convl_gu*f1;   bnorm_sq = bnorm_sq+ sum(f1.^2)*dx; 
         gu_fN2(:,nn)   = convl_gu(:,1:2:end)*f1(1:2:end); 
    end 

    rho_val  = mean(rhoN,2);    
    indr = find(rho_val>0);  rho_val = rho_val(indr);  % remove those zero weight points 
    r_seq   = r_seq(indr);    
    regressionData.rho_val = rho_val/sum(rho_val)/dx;   % exploration measure  
    Gbar     = mean(GbarN,3)*dx;     regressionData.Gbar = Gbar(indr,indr);      % integral kernel in the operator L_G       *dr*ds
    gu_f     = mean(gu_fN,2)*dx;     regressionData.gu_f = gu_f(indr);      % b                                      *dr
    gu_f2    = mean(gu_fN2,2)*dx;    regressionData.gu_f2= gu_f2(indr);     % b2 
	regressionData.bnorm_sq   = bnorm_sq/N;
	regressionData.bdry_width = bdry_width; 
	regressionData.r_seq      = r_seq;  
	regressionData.data_str   = data_str;

end  


end
