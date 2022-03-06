function  [kernelInfo, obsInfo]= load_settings(N, u_Type, u_supp, dx, kernel_type, example_type)

xDim   = 1;       % dimension

%% set example types; define the function in convolution with phi
% L_phi[u](x) = \int_y phi(|y|)*fun_g[u](x,y) dy
%             = \sum_r phi(r)* [ fun_g[u](x,x+r) + fun_g[u](x,x-r) ] dr

%if ~exist('example_type','var'); example_type = 'nonlocal'; end

% example_type = 'LinearIntOpt';  % {'LinearIntOpt','nonlocal','nonlinearOpt'};
temp = strsplit(example_type, '_');
switch strtok(example_type, '_')
    % 1. g[u](x,y) =  u(x+y)          integral kernel: 1D, Laplace Transform
    case 'LinearIntOpt'
        fun_g = @(u,Du,x,y) u(x+y);
        fun_g_vec = @(u,Du,xInd,rInd) u(xInd+rInd);
        % 2. g[u](x,y) =  [u(x+y)-u(x)]   nonlocal kernel
    case 'nonlocal'
        fun_g = @(u,Du,x,y) u(x+y)-u(x);
        fun_g_vec = @(u,Du,xInd,rInd) u(xInd+rInd)- u(xInd);
        % 3. g[u](x,y) = div(u(x+y))u(x)  particle interaction
    case 'nonlinearOpt'
        fun_g = @(u,Du,x,y)  Du(x+y).*u(x);
        fun_g_vec = @(u,Du,xInd,rInd) Du(xInd+rInd).*u(xInd);
    case 'mfOpt'
        fun_g = @(u,Du,x,y)  Du(x+y).*u(x) + u(x+y).*Du(x);
        fun_g_vec = @(u,Du,xInd,rInd) Du(xInd+rInd).*u(xInd)+ u(xInd+rInd).*Du(xInd);
    case 'ParticleSystem'
        fun_g = @(u,Du,x,y) u(x+y)-u(x);
        fun_g_vec = @(u,Du,xInd,rInd) u(xInd+rInd)- u(xInd);
        particles.d = str2double(temp{2}(2));
        particles.pN = str2double(temp{3}(2:end));
        particles.mu = mean(u_supp);
        particles.sigma = (u_supp(2) - particles.mu)/4;
        obsInfo.particles = particles;
    case 'MeanField'
        fun_g = @(x) 0*x;
        fun_g_vec = @(x) 0*x;
    otherwise; fun_g = @(u,Du,x,y)  ones(size(x)); fun_g_vec = @(u,Du,xInd,rInd) ones(size(xInd));
end
obsInfo.example_type = example_type;
obsInfo.fun_g = fun_g;
obsInfo.fun_g_vec = fun_g_vec;

kernelInfo.d = xDim;
%% true kernel K: R+ -->R
if ~exist('kernel_type','var'); kernel_type = 'Gaussian'; end
% kernel_type = 'Gaussian';                       % sinkx,  GaussianPN, FracLap
switch kernel_type
    case 'sinkx';    k = 4;
        K_true = @(r) sin(k*r).*(r>=0).*(r<=3);
        threshold = 1e-8;   % threshold for support
        kernel_str = 'sinkx';
    case 'Gaussian';    s  = 0.75; mu = 3;
        K_true = @(r) exp(-0.5/s/s * (r - mu).^2)/sqrt(2*pi)/s;
        kernel_str = [kernel_type,'_mean_',num2str(mu), '_std_', num2str(s)];
        threshold = 1e-10;
    case 'GaussianPN';  s1 = 1; mu1 = 5; s2= 2; mu2 = 0;
        K1_true = @(r) exp(-0.5/s1/s1 * (r - mu1).^2)/sqrt(2*pi)/s1;
        K2_true = @(r) exp(-0.5/s2/s2 * (r - mu2).^2)/sqrt(2*pi)/s2;
        K_true  = @(r) K1_true(r) - K2_true(r);
        kernel_str = [kernel_type,'_mean1_',num2str(mu1), '_std1_', num2str(s1),'_mean2_',num2str(mu2), '_std2_', num2str(s2)];
        threshold = 1e-8;
    case 'FracLap';     s = 0.5; d = xDim; c_ds = (4^s*gamma(d/2+s))./((pi^(d/2)) * abs(gamma(-s)));
        K_true = @(r) c_ds * ((1./r.^(d+2*s)) .* (r>0.1).*(r<=3) + (1./0.1.^(d+2*s)) .* (r<=0.1))+0*0.00842*(r>3);
        kernel_str = [kernel_type,'_d_',num2str(d), '_s_', num2str(s)];
        kernel_str = strrep(kernel_str, '.', '');
        kernelInfo.s = s;  kernelInfo.c_ds= c_ds;
        threshold = 1e-4;
end
kernelInfo.K_true = K_true;
kernelInfo.kernel_type = kernel_type;
kernelInfo.kernel_str = kernel_str;

obsInfo.threshold = threshold;
obsInfo.u_supp    = u_supp;


%% function handle of u and f, in terms of basis functions
% known u_i, f_i will be evaluated at points x_j later
if strcmp(example_type, 'classicalReg') == 1
    obsInfo.u{1} = @(x) x;
    obsInfo.div_u{1} = @(x) ones(size(x));
    obsInfo.u_str = '';
else
    switch u_Type
        case 'Bspline'
            degree   = 1; % 1, 2, 3...
            knotsAug = 0;
            knots    = linspace(u_supp(1), u_supp(end), N+degree)';
            [v, vknots] = Basis_bspline(degree, knots, knotsAug);
            obsInfo.u = v;
            obsInfo.u_knots = vknots;obsInfo.u_degree = degree;
            obsInfo.u_str   = [u_Type,'_',num2str(knots(1)), '_',num2str(knots(end)), '_N_',num2str(numel(v)),'_deg_',num2str(degree), '_knotsAug_',num2str(knotsAug)];
        case 'Fourier'
            odd_even = 'sine';    % 'sine', 'cosine'  % use sine or cosine series, or use both;
            obsInfo.u_str   = [u_Type,'_N_',num2str(N), odd_even];
            [v,basisInfo] = Basis_fourier(N, odd_even, u_supp);
            obsInfo.u = v;
            obsInfo.div_u = basisInfo.dbasis_funs;
            obsInfo.laplacian_u = basisInfo.ddbasis_funs;
        case 'Particles_Gaussian'
            odd_even = 'sincos';    % 'sine', 'cosine'  % use sine or cosine series, or use both;
            obsInfo.u_str   = [u_Type,'_N_',num2str(N), odd_even];
            [v,~] = Basis_fourier(N, odd_even, u_supp);
            obsInfo.u = v;
            obsInfo.div_u = v;
        case 'GaussianMixture'
            temp = cell(N, 1);
            lb = u_supp(1);
            ub = u_supp(2);
            md = (lb + ub)/2;
            wd = md-lb;
            lb_mu = md - wd/2;
            ub_mu = md + wd/2;
            std_max = wd/6;
            mixture_num = 3;
            for i = 1:N
                mu = unifrnd(lb_mu, ub_mu, [1, mixture_num]);
                std = unifrnd(0, std_max, [1, mixture_num]);
                
                density = @(x) 0*x;
                for j = 1:mixture_num
                    density = @(x) density(x) + normpdf(x, mu(j), std(j));
                end
                temp{i} = density;
            end
            obsInfo.u = temp;
            obsInfo.u_str = 'GaussianMixture';
    end
end
%% add u(x) = x for case nonlinear opt
if strcmp(u_Type, 'Fourier')
    if strcmp(example_type, 'nonlinearOpt') == 1 || strcmp(example_type, 'mfOpt') == 1
        obsInfo.u{N+1} = @(x) x.*(x>=u_supp(1)).*(x<=u_supp(2));
        obsInfo.div_u{N+1} = @(x) ones(size(x)).*(x>=u_supp(1)).*(x<=u_supp(2));
        obsInfo.laplacian_u{N+1} = @(x) zeros(size(x));
        obsInfo.u_str   = strrep(obsInfo.u_str,['N_',num2str(N)],['N_',num2str(N+1)]); % sttrep(str,old,new) % 
    end
end
% default x_mesh
range  = 40;    % the support of the kernel should be in [0,range]; will be refined below to be adaptive to kernel
if strcmp(example_type, 'classicalReg') == 1 || strcmp(obsInfo.example_type, 'nonlinearOpt') == 1 ||  strcmp(obsInfo.example_type, 'mfOpt') == 1
    xmin   = u_supp(1);
    xmax   = u_supp(end);
else
    xmin   = u_supp(1) - range;
    xmax   = u_supp(end) + range;
end
x_mesh = [xmin:dx:xmax];
x_mesh_str = [num2str(xmin), '_', num2str(dx),'_', num2str(xmax)];
x_mesh_str = strrep(x_mesh_str, '.','_');
obsInfo.x_mesh     = x_mesh;
obsInfo.x_mesh_dx  = dx;
obsInfo.x_mesh_str = x_mesh_str;

% function f with true kernel: in analytically format
f = cell(size(obsInfo.u)); intgrand = cell(size(f));

if ~(strcmp(example_type, 'MeanField')) && ~(strcmp(example_type, 'Homogenized'))
    for  i = 1:numel(f)
        if strcmp(example_type, 'classicalReg') == 1
            f{1} = @(x) K_true(x);
        elseif strcmp(example_type, 'nonlinearOpt') == 1 || strcmp(example_type, 'mfOpt') ==1
            intgrand{i} = @(x,y) K_true(abs(y)).* fun_g(obsInfo.u{i},obsInfo.div_u{i},x,y); 
            f{i} =@(x) obsInfo.laplacian_u{i}(x) + quadgk(@(y)intgrand{i}(x,y), -Inf, Inf,'Waypoints', obsInfo.x_mesh, 'MaxIntervalCount', 10^10);
        else
            intgrand{i} = @(x,y) K_true(abs(y)).* fun_g(obsInfo.u{i},obsInfo.div_u{i},x,y);
            %f{i} = @(x) integral(@(y)intgrand(x,y), xmin, xmax);   % to explain this algorithm in note: global adaptive quadrature
            % why it does not work --- the f should be smooth.
            f{i} =@(x) quadgk(@(y)intgrand{i}(x,y), -Inf, Inf,'Waypoints', obsInfo.x_mesh, 'MaxIntervalCount', 10^10);  
        end
    end
end
obsInfo.f = f;
%obsInfo.intgrand = intgrand;

%obsInfo = conv_kernel(u_Type,N,supp_u,K_true,dx)

fprintf('Example type: %s \n Kernel type: %s\n ',example_type,kernel_type);

obsInfo = Effective_kernel_adaptive_x_mesh(obsInfo);

end


function obsInfo = Effective_kernel_adaptive_x_mesh(obsInfo)
% kernel-adaptive mesh: Use the true kernel in ananlytical form to detect the actual range of
% interaction, so as avoid no-used x_mesh. This can significantly save computational cost.
u_supp = obsInfo.u_supp;
dx        = obsInfo.x_mesh_dx;

switch obsInfo.example_type
    case {'nonlinearOpt','mfOpt'}
        supp_H     = [0 u_supp(end) - u_supp(1)+dx];
        % Set kernel adaptive x_mesh:
        range = supp_H(2);
        xmin   = obsInfo.u_supp(1) - range;
        xmax   = obsInfo.u_supp(end) + range;
        x_mesh = xmin:dx:xmax;
        x_mesh_str = [num2str(xmin), '_', num2str(dx),'_', num2str(xmax)];
        x_mesh_str = strrep(x_mesh_str, '.','_');
        dx         = x_mesh(2) - x_mesh(1);
        obsInfo.x_mesh     = x_mesh;
        obsInfo.x_mesh_str = x_mesh_str;
        obsInfo.x_mesh_dx  = dx;
    case {'classicalReg', 'MeanField', 'Homogenized'}
        supp_H     = [u_supp(1)-dx u_supp(end)+dx];
    otherwise
        threshold = obsInfo.threshold;
        u         = obsInfo.u{1};
        f         = obsInfo.f{1};   % it uses the true kernel in analytical form
        default_x_mesh = obsInfo.x_mesh;
        % find xmin, xmax such that f(xmin), f(xmax) < threshold
        fval = arrayfun(f,default_x_mesh); % figure; plot(fval)
        ind_min = find(abs(fval)>=threshold, 1, 'first') - 1;
        ind_max = find(abs(fval)>=threshold, 1, 'last') + 1;
        ind_min = max(ind_min, 1);
        ind_max = min(ind_max, numel(fval));
        % range of K from data u, f
        R = max([default_x_mesh(ind_max) - u_supp(2), u_supp(1) - default_x_mesh(ind_min)]);
        supp_H = [0 1.1*R];
        % Set kernel adaptive x_mesh:
        range = supp_H(2);
        xmin   = obsInfo.u_supp(1) - 2*range;
        xmax   = obsInfo.u_supp(end) + 2*range;
        x_mesh = xmin:dx:xmax;
        x_mesh_str = [num2str(xmin), '_', num2str(dx),'_', num2str(xmax)];
        x_mesh_str = strrep(x_mesh_str, '.','_');
        dx         = x_mesh(2) - x_mesh(1);
        obsInfo.x_mesh     = x_mesh;
        obsInfo.x_mesh_str = x_mesh_str;
        obsInfo.x_mesh_dx  = dx;
end

obsInfo.supp_H     = supp_H;
obsInfo.delta      = round(supp_H(end),3);
fprintf('\n Kernel support range [%2.2f, %2.2f]\n',supp_H(1), supp_H(2));

end

