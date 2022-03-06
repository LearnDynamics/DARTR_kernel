function [x_reg,lambda_opt] = Lcurve_with_Norm_pinv(A, b, B, titl, plotON)
% L-curve regularization with a given norm and specified range for parameter

%% old
% N = 1000;
% lambda_seq = 10.^(linspace(-30, 10, N));
% len = length(lambda_seq);
% E = zeros(len,1);
% R = zeros(len,1);
% const = b'*(pinv(A) * b);
% for ll = 1:len
%     lambda1 = lambda_seq(ll);
%     x_l     = pinv(A+lambda1*B) * b;
%     E(ll) = x_l'*A*x_l -2*b'*x_l + const;
%     R(ll) = x_l'*B*x_l;
% end

%% new

% Some observations:
% 1. c'Ac - 2b'c + b'A^-1b = (Ac - b)'A^-1(Ac-b)
% 2. (A + lambda B)c = b  ====> Ac-b = -lambda * B * c
N = 1000;
lambda_seq = 10.^(linspace(-16, 10, N));
len = length(lambda_seq);
E = zeros(len,1);
R = zeros(len,1);

A_inv = pinv(A);
A2 = sqrtm(A);
b2 = pinv(A2)*b;
for ll = 1:len
    lambda1 = lambda_seq(ll);
    x_l     = pinv(A+lambda1*B)*b;% lsqminnorm(A+lambda1*B, b);
%     E(ll) = (A*x_l-b)'*lsqminnorm(A, A*x_l-b);
%     E(ll) = x_l'*A*x_l -2*b'*x_l + b'*(pinv(A) * b);
    E(ll) = norm(A2*x_l - b2); %This is the best method
    R(ll) = x_l'*B*x_l;
end


% from eig(Abar),eig(Abar,Bmat)
% N = 200;
% lambda_seq = lambda_seq_from_regMats(A_bar,Bmat,N); 
%A_bar_sq = A_bar*A_bar;
%A_bar_b_bar = A_bar*b_bar;

if 0
%%
    figure;
    subplot(2,1,1);scatter(log10(lambda_seq), log10(R));xlabel log_{10}(\lambda);ylabel log_{10}(R);
	subplot(2,1,2);scatter(log10(lambda_seq), log10(E));xlabel log_{10}(\lambda);ylabel log_{10}(E);
end


%% L-curve
if ~exist('curvatureType','var');  curvatureType = 'fit_curve'; end
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