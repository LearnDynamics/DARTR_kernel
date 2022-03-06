%%  learn nonlocal kernel from data 
%{
Given data (u_i,f_i), we fit a nonlocal kernel to it: 
                L_phi[u_i] = f_i;     
where  L_phi[u](x) = \int phi(|y|)*g[u](x,x+y) dy 
                   = \sum_r phi(r)* [ g[u](x,x+r)+ g[u](x,x-r) ] dr 
1. g[u](x,y) =  u(x+y)          integral kernel: 1D, Laplace Transform 
2. g[u](x,y) =  [u(x+y)-u(x)]   nonlocal kernel
3. g[u](x,y) = div(u(x+y))u(x)  mfOpt 
% Here we consider only 1D; 
%}
clear all; close all; restoredefaultpath;  
add_mypaths;
rng(1);
plotOn = 1; saveON = 0; 

%% 0 setttings
dx       = 0.05;                          % space mesh size
N        = 2;                             %  number of data pairs (u_i,f_i)
u_Type      = 'Fourier';   supp_u = [-1*pi 1*pi];  % data u
kernel_type = 'sinkx';       % Gaussian, sinkx, FracLap
example_type = 'mfOpt';      % {'LinearIntOpt','nonlocal','mfOpt'};


SAVE_DIR     = [SAVE_DIR,example_type,'/']; if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end 

[kernelInfo, obsInfo]  = load_settings(N, u_Type, supp_u, dx,kernel_type,example_type);
noise_ratio            = 1;  obsInfo.noise_ratio  = noise_ratio;                                 % noise to signal ratio

run_1setting;
