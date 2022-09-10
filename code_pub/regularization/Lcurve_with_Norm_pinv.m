function [x_reg,lambda_opt] = Lcurve_with_Norm_pinv(A, b, B, titl, plotON)
% L-curve regularization with a given norm and specified range for parameter

%% evaluate error/loss and regularization term
          % lambda in range adapative to eigenvalues
          % loss evaluated by norm(sqrt(A)x-b): does not need to estimate the constant 
          % x_lambda solved by pinv or lsqmininorm
[~, eigL]  = eig(A,B);  % generalized eigenvalue   A V = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
eigL       = real(diag(eigL));

% Some observations:
% 1. c'Ac - 2b'c + b'A^-1b = (Ac - b)'A^-1(Ac-b)
% 2. (A + lambda B)c = b  ====> Ac-b = -lambda * B * c
N = 1000;
if plotON==1; N=200; end
lambda_seq = 10.^(linspace(log10(min(eigL)+1e-16)-3, log10(max(eigL)), N));
%lambda_seq = 10.^(linspace(-16, 10, N));
len = length(lambda_seq);
E = zeros(len,1);    % error/loss 
R = zeros(len,1);    % regularization norm^2

A2 = sqrtm(A);   b2 = pinv(A2)*b;    % so that E = norm(A2*x_l - b2);

compare_computE(lambda_seq,A,B,b); 

for ll = 1:len
    lambda1 = lambda_seq(ll);
    x_l     = pinv(A+lambda1*B)*b; % lsqminnorm(A+lambda1*B, b);
%     E(ll) = (A*x_l-b)'*lsqminnorm(A, A*x_l-b);
%     E(ll) = x_l'*A*x_l -2*b'*x_l + b'*(pinv(A) * b);
    E(ll) = norm(A2*x_l - b2); % This is the best method
    R(ll) = x_l'*B*x_l;
end
R = sqrt(R); 


if 0   % plot the scattered L-curve points
    figure;
    subplot(2,1,1);scatter(log10(lambda_seq), log10(R));xlabel log_{10}(\lambda);ylabel log_{10}(R);
	subplot(2,1,2);scatter(log10(lambda_seq), log10(E));xlabel log_{10}(\lambda);ylabel log_{10}(E);
end


%% L-curve
if ~exist('curvatureType','var');  curvatureType = '3ptCircle'; end
switch curvatureType
    case 'fit_curve'
        lambda_opt = Opt_lambda_Quanjun(E, R, lambda_seq, plotON, titl);
    case '3ptCircle'      % use 3point circle radius
        lambda_opt = Opt_lambda_curvature3pt(E,R, lambda_seq, plotON,titl);
    otherwise
        lambda_opt = opt_lambda_curvature(E,R, lambda_seq, plotON);
end
x_reg =  real(pinv(A+lambda_opt*B)*b); % (A+lambda_opt*B)\b;
end


function compare_computE(lambda_seq,A,B,b)
% compare two ways computing E: by c'*A*c - 2*c'*b + const. 
% conclusion: they seem similar, but could have big effect when lambda is small.  
len = length(lambda_seq);
E2 = zeros(len,1);    % error/loss 
R2 = zeros(len,1);    % regularization norm^2
E  = zeros(len,1);
% Approach 2: loss evaluated by norm(sqrt(A)x-b): does not need to estimate the constant 
% Approach 1: direct evaluation of loss: with const estimated 
A2    = sqrtm(A); 
b2    = pinv(A2)*b;
const = b'*(pinv(A) * b);  
for ll = 1:len
    lambda1 = lambda_seq(ll);
    x_l     = pinv(A+lambda1*B)*b; % lsqminnorm(A+lambda1*B, b);
%     E(ll) = (A*x_l-b)'*lsqminnorm(A, A*x_l-b);
%     E(ll) = x_l'*A*x_l -2*b'*x_l + b'*(pinv(A) * b);
    E2(ll) = norm(A2*x_l - b2)^2; % This is the best method
    E(ll)  = x_l'*A*x_l -2*b'*x_l + const;
    R2(ll) = x_l'*B*x_l;
end
figure; 
loglog(lambda_seq,E,lambda_seq,E2,'-.','linewidth',1);xlabel('lambda'); ylabel('loss values');
legend('E with estimated const','E= norm(sqrt(A)-b)'); 
end