function  [obsInfo,ux_val,fx_val] = generateData(kernelInfo, obsInfo, SAVE_DIR,saveON)
% generate data without noise

if ~exist('saveON','var'); saveON = 0; end
filename = [SAVE_DIR, 'Data_', kernelInfo.kernel_str, '_ui_', obsInfo.u_str,'_xj_', obsInfo.x_mesh_str, '.mat'];


if ~exist(filename, 'file')
    f      = obsInfo.f;
    u      = obsInfo.u;
    x_mesh = obsInfo.x_mesh;
    J      = numel(x_mesh)^kernelInfo.d; N = numel(u);
    % dist_x = abs(x_mesh' - x_mesh); % pairwise dist: dist_x(j1, j2) = |x_j1 - x_j2|
    ux_val = zeros(N, J); % ux_val(n,j) = u_n(x_j)
    % dux_val = zeros(N, J, J); % dux_val(n,j1,j2) = u_n(x_j1) - u_n(x_j2)
    fx_val = zeros(N, J); % fx_val(n,j) = f_n(x_j) --> b_mat
    %{
    if kernelInfo.d == 2
        [meshx1, meshx2]= meshgrid(x_mesh, x_mesh);
        for n = 1:N
        ux_val(n,:) = reshape(arrayfun(u{n}, meshx1, meshx2), 1,[]);
        dux_val(n,:,:) = ux_val(n,:)' - ux_val(n,:);
        fx_val(:,n) = reshape(arrayfun(f{n}, meshx1, meshx2), 1,[]);
        end
    end
    %}
    switch obsInfo.example_type
        case {'LinearIntOpt','nonlocal','nonlinearOpt','mfOpt'}
            for n = 1:N
                ux_val(n,:) = arrayfun(u{n}, x_mesh);
                temp        = arrayfun(f{n}, x_mesh);
                fx_val(n,:) = temp';
            end
        case {'MeanField', 'Homogenized'}
            for n = 1:N
                ux_val(n,:) = arrayfun(u{n}, x_mesh);
                fx_val(n,:) = R_phi_u(u{n}, kernelInfo.K_true, x_mesh, obsInfo.example_type);
            end
        otherwise
            if contains(obsInfo.u_str, 'Particles')
                p = obsInfo.particles;
                d = p.d;
                pN = p.pN;
                mu = p.mu;
                sigma = p.sigma;
                phi = kernelInfo.K_true;
                ux_val = randn(N, pN, d)*sigma+mu;
                fx_val = zeros(N, pN, d);
                fx_val = particles_f(phi, ux_val);
            end    
    end
    
    
    % obsInfo.meshed_data.dist_x = dist_x;
    %    obsInfo.meshed_data.ux_val = ux_val;
    %obsInfo.meshed_data.dux_val = dux_val;
    %    obsInfo.meshed_data.fx_val = fx_val;
    
    %% noise range to be reasonable: depend on mesh, depend on f
    % noiseL2^2 = noise_std^2 *sum(length(x_mesh))*dx
    % fL2^2     = sum(f(:,1).^2)*dx
    
    
    dx = obsInfo.x_mesh_dx;
    f_L2norm2   = sum(fx_val.^2,2)*dx;
    f_mean      = mean(fx_val,2);
    f_std       = std(fx_val,0,2);
    noise_std_upperBd = min( [ min(f_mean+f_std), sqrt(min(f_L2norm2)/sum(length(x_mesh))*dx)]);
    obsInfo.noise_std_upperBd = noise_std_upperBd;
    obsInfo.f_L2norm  = f_L2norm2;
    obsInfo.f_mean    = f_mean;
    obsInfo.f_std     = f_std;
    
    if contains(obsInfo.u_str, 'Particles')
        obsInfo.noise_std_upperBd = mean(squeeze(noise_std_upperBd), 2);
        obsInfo.f_L2norm  = mean(squeeze(f_L2norm2), 2);
        obsInfo.f_mean    = mean(squeeze(f_mean), 2);
        obsInfo.f_std     = mean(squeeze(f_std), 2);
        obsInfo.noise_std_upperBd = min( [ min(obsInfo.f_mean+obsInfo.f_std), sqrt(min(obsInfo.f_L2norm)/sum(length(x_mesh))*dx)]);
    end
    
    if saveON==1;     save(filename,  'obsInfo','ux_val','fx_val'); end
else
    load(filename,  'obsInfo','ux_val','fx_val'); %
end

end