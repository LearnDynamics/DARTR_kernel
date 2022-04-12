%% Demonstration of the regularization
% given Q, aplpha, to estmate f 
%{
  Q(t)     = sum_{i=1}^nv alpha(t) f(vi) vi* exp(-vi*t) 
  alpha(t) = sumsum(J), with J(t) randomly sampled. 
%}

folder = fileparts(which(mfilename));
addpath(genpath(folder));

v_min = 1/500; 
v_max = 0.6; 
vgrid = linspace(v_min,v_max,100); 
dv    = vgrid(2)- vgrid(1); 

ftrue = @(v) exp(sin(10*v)); 
gvt   = @(v,t) v.* exp(-v*t); 

nv    = length(vgrid); 
m     = 1000;  dt = 0.1;  
tseq  = dt*(1:m); 
J     = rand(1,m); 
alphat= cumsum(J); 
fseq  = ftrue(vgrid); 

%% Data of Qseq
Qseq  = zeros(1,m); 
for t=1:m
    gvt_vseq = gvt(vgrid,tseq(t)); 
    Qseq(t)  = alphat(t)*sum(fseq.*gvt_vseq);
end


%% regression to estimator f
% get A, b  exploration measure rho
A = zeros(nv,nv);  b= zeros(nv,1); rhof = zeros(nv,1);
for i=1:nv
    for j=i:nv
        vivj  = vgrid(i) * vgrid(j); 
        viPvj = vgrid(i) + vgrid(j); 
        A(i,j) = vivj*sum(alphat.^2.*exp(-viPvj*tseq))*dt;
        A(j,i) = A(i,j);
    end
    vi   = vgrid(i); 
    b(i) = vi*sum(alphat.*exp(-vi*tseq).*Qseq)*dt; 
    rhof(i) = vi*sum(alphat.*exp(-vi*tseq))*dt;  
end
rhof = rhof/sum(rhof);


Bmat = diag(rhof); titl ='RKHS regu'; plotON =1; 
[x_reg,lambda_opt] = Lcurve_sidaRKHS_lsq2(A,b,Bmat,titl,plotON);

x_no_reg = A\b; 

figure; 
plot(vgrid,fseq,'k:','linewidth',2); hold on;
plot(vgrid,x_reg,'b--','linewidth',1); 
plot(vgrid,x_no_reg,'r-.','linewidth',1); 
xlabel('v'); ylabel('f(v)'); 
legend('True','Estimated DARTR','Estimated A\\b')
ylim([min(fseq)*0.9,max(fseq)*1.2])

ftsz = 14; %mksz = 12; % if save to -dpdf
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
set_positionFontsAll;
print('Demo_deconvolution.pdf','-dpdf', '-bestfit');






