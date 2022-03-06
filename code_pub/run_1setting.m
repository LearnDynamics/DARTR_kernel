% function [Est_vector,Est_ftn_dimSeq] = run_1setting(kernelInfo,obsInfo,SAVE_DIR,saveON,plotOn)

%% 1. Generate data on x_mesh (adaptive to u, f) and pre-process data
[obsInfo,ux_val,fx_val]   = generateData(kernelInfo, obsInfo, SAVE_DIR,saveON);    % get observation data
noise_std_upperBd = max(round(obsInfo.noise_std_upperBd,2),0.001);
nsr               = obsInfo.noise_ratio;
noise_std         = nsr *noise_std_upperBd;  obsInfo.noise_std = noise_std;% noise added to fx_val
fx_val            = fx_val + noise_std.*randn(size(fx_val));

dx        = obsInfo.x_mesh_dx; 
data_str  = [obsInfo.example_type,kernelInfo.kernel_type,obsInfo.u_str,obsInfo.x_mesh_str,sprintf('NSR%1.1f_dx%1.4f_',nsr,dx)];

if strcmp(obsInfo.example_type, 'classicalReg') == 1 %|| strcmp(obsInfo.example_type, 'nonlinearOpt') == 1
    bdry_width = 0; 
    r_seq      = obsInfo.x_mesh;    
else
    bdry_width = floor((obsInfo.delta+1e-10)/dx);  % boundary space for inetration range with kernel 
    r_seq      = dx*(1:bdry_width);
end
Index_xi_inUse = (bdry_width+1): (length(obsInfo.x_mesh) - bdry_width); % index that x_i in use for L_phi[u]


%% 2 pre-process data, get all data elements for regression:   ****** a key element to significantly reduce computational cost
div_u = zeros(size(ux_val));
if ~contains(obsInfo.example_type, 'ParticleSystem')
    for n = 1:length(ux_val(:,1))
        div_u(n,:) = obsInfo.div_u{n}(obsInfo.x_mesh);
    end
end
normalizeOn = ~strcmp(obsInfo.example_type, 'classicalReg');
fun_g_vec = obsInfo.fun_g_vec;
regressionData = getData4regression(ux_val,div_u,fx_val,dx,obsInfo,bdry_width,Index_xi_inUse,r_seq,data_str,normalizeOn);
clear ux_val fx_val;

N_r_seq        = numel(regressionData.r_seq);
%% 3. Estimator: vector estimator and function estimator
Bmat_type  = 'L2rho';    % 'L2rho', 'Unif_rho','Lebsegue': compute Bmat = <phi_i, phi_j> using measure rho, Lebesgue or discrete-uniform
plot4paper = 1; 

normType = {'l2','L2','RKHS'};   % ,'H1' to compare H1, otherRKHS?
%% vector estimator: value of the kernel on meshes
 [Est_vector] = vector_kernel_run2(kernelInfo, regressionData, SAVE_DIR, plotOn,saveON,Bmat_type,normType,plot4paper);

%% 3. function estimator: with basis functions. 

