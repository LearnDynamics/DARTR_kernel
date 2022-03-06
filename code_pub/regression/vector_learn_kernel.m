function [Est,h_fig_coef] = vector_learn_kernel(regressionData, plotOn,Bmat_type,normType)
% estimate kernel vector with regularization; coefficient c gives the values of kernel 
% test four regularization strategies: truncated SVD with L2 or RKHS; L-curve with L2 or RKHS; 
% 
%% 1 compute matrices for regression, and rho and Gbar
rho_val = regressionData.rho_val; 
ind_rho = find(rho_val>0);  rho_val = rho_val(ind_rho);
r_seq   = regressionData.r_seq(ind_rho);         % when r_seq is non-uniform, use dr = r_seq(2:end) - r_seq(1:end-1).       
dr      = r_seq(2)-r_seq(1); 
Abar    = regressionData.Gbar(ind_rho,ind_rho)*dr*dr;     % A(i,j) = <L_{phi_i}[u], L_{phi_j}[u]> = sum_{r,s} phi_i(r)phi_j(s) G(r,s) drds =>> A(i,j) = G(r_i,r_j) *dr*dr
bbar    = regressionData.gu_f(ind_rho)*dr;        % b(i)   = <L_{phi_i}[u], f> = int_{r} phi_i(r) sum_x convl_gu(r,x)f(x)dx dr     =>> b(i)   = convl_gu_f(r_i) *dr
bbar2   = regressionData.gu_f2(ind_rho)*dr; 
B_L2rho = diag(rho_val())*dr;               % B_L2rho(i,j) = <phi_i,phi_j>_rho = diag(rho_val)*dr
BmatUnif= eye(length(bbar))/(dr*length(r_seq));              % Bmat         = <phi_i,phi_j>_Lebesgue = eye/(dr*length(r_seq) )
b_norm_sq = regressionData.bnorm_sq; 

%% 2. estimators: coefficient c gives the values of kernel
% regularization:  SVD with l2, L2, and RKHS; Lcurve with l2, L2 and RKHS
% When B = eye:      l2 = L2 =RKHS;  
% when B = B_L2rho:  l2 = above, L2 = L2(rho), RKHS = range(Abar);   
switch Bmat_type     % 'L2rho', 'Unif_rho','Lebsegue'
    case 'L2rho';     Bmat = B_L2rho;      
    case 'Unif_rho';  Bmat = BmatUnif;        % Bmat computed with uniform measure on discrete set 
    case 'Lebsegue';  Bmat = BmatLebesgue;     % Bmat computed with Lebesgue measure 
          fprintf('Vector estimator: BmatBmatLebesgue is used, but should not apply\n'); 
end
tol = 1e-20; 
Bmat   = Lift_smallestEigen(Bmat,tol);  % replace the small eigval by tol; to avoid numerical issues in the generalized eigenvalue problem
%% Compare 3 RKHS and H1
if ~exist('normType','var')
    normType = {'l2','L2','RKHS'};   % ,'H1' to compare H1, otherRKHS?
end
[~,~,cLcurveAll] = reg_Lcurve_3in1(Abar,bbar,Bmat,plotOn, normType);
cLcurveAll.normType = normType;     
%[c_plain,cLcurveAll] = regu_Lcurve_3RKHS_H1(Abar,bbar,Bmat,plotOn);
Est.matInfo.Abar  = Abar;  
Est.matInfo.bbar  = bbar;              Est.matInfo.bbar2   = bbar2; 
Est.matInfo.Bmat  = Bmat;    
Est.matInfo.b_norm_sq = b_norm_sq;  
Est.rho_val       = rho_val; 
Est.r_seq         = r_seq;
Est.cLcurveAll     = cLcurveAll;   

h_fig_coef =[]; 
% plot estimated kernel
if plotOn == 1
    h_fig_coef = figure; xi = r_seq; c_plain = Abar\bbar; 
    plot(xi, c_plain, 'c','Displayname', 'LSE');hold on;
    Extreme_val = zeros(length(normType),2);  % find the max and min of estimators
    for i=1:length(normType)
        fieldname = ['creg_',normType{i}]; estimator = cLcurveAll.(fieldname);
        plot(xi, estimator, '-.','Displayname', normType{i}, 'Linewidth', 2);
        Extreme_val(i,:) = [max(estimator), min(estimator)];
    end
    ylim([ 0.8*min(Extreme_val(:,2)),1.2*max(Extreme_val(:,1))]);
    % % super-impose rho
    normalizer = max(estimator)/2/max(rho_val); rho_length = length(rho_val);
    nx = length(xi); gap = ceil(nx/rho_length); indx = 1:gap:nx; xirho = xi(indx); % sparser index
    area(xirho,rho_val(indx)*normalizer, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeAlpha', 0,'Displayname', '\rho');
    hold off; legend();title('Coefficient Estimators'); 
    
   figure(101); clf; subplot(121);  plot(xi,rho_val,'b-.','linewidth',2);
        xlabel('r'); ylabel('Exploration density'); % ytickformat('%1.4f') % xtickformat('$%,.0f')
        Gbar_rho =  Abar./(rho_val*rho_val');
        [~,eigG] = eig(Abar); [~,eigGrho] = eig(Gbar_rho,diag(rho_val));
        eigG = sort(diag(eigG),'descend'); eigGrho = sort( diag(eigGrho),'descend');
        subplot(122);
        semilogy(eigG,'b-.','linewidth',1); hold on;
        semilogy(eigGrho,'r-','linewidth',1);
        legend('Eigenvalue Gbar','Eigenvalue Gbar-rho L2(rho)');
end


end


    
    
