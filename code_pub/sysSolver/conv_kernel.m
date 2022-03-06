function obsInfo = conv_kernel(u_Type,N,supp_u,K_true,dx)
% % function handle of u and f, in terms of basis functions  
% Q: Is it faster than direct solver without using such a handle? 

%% known u_i, f_i at points x_j 
% u_i
%N = 3; 
%u_Type = 'Bspline';  % Fourier or Bspline
%supp_u = [0 20];
switch u_Type
    case {'Bspline'}
        degree   = 3; % 1, 2, 3...
        knotsAug = 0;
        knots = linspace(supp_u(1), supp_u(end), N+degree)';
        [u, uInfo] = Basis_bspline(degree,knots,knotsAug);
        obsInfo.u = u;
        obsInfo.u_knots = uInfo.knots;
        obsInfo.u_str   = [u_Type,'_',num2str(knots(1)), '_',num2str(knots(end)), '_N_',num2str(numel(u)),'_deg_',num2str(degree), '_knotsAug_',num2str(knotsAug)];
    case 'Fourier'
        odd_even = 'sincos';    % 'sine', 'cosine'  % use sine or cosine series, or use both;  
        obsInfo.u_str   = [u_Type,'_N_',num2str(N), odd_even]; 
        [u,~] = Basis_fourier(N, odd_even);
        obsInfo.u = u;
end

% f_i
f = cell(size(u));
for  i = 1:numel(u)
    intgrand = @(x,y) K_true(abs(y-x)).*(obsInfo.u{i}(y) - obsInfo.u{i}(x));
    f{i} = @(x) integral(@(y)intgrand(x,y), -Inf, Inf);
end

obsInfo.f = f;

%% Set x_mesh 
xmin = obsInfo.u_knots(1) - supp_H(end); xmax = obsInfo.u_knots(end) + supp_H(end);

x_mesh = [xmin:dx:xmax]; 
x_mesh_str = [num2str(xmin), '_', num2str(dx),'_', num2str(xmax)];
x_mesh_str = strrep(x_mesh_str, '.','_');
obsInfo.x_mesh = x_mesh;
obsInfo.x_mesh_dx = dx;
obsInfo.x_mesh_str = x_mesh_str;
