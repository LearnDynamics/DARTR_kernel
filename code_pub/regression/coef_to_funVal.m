function [K_est_val,structName] = coef_to_funVal(coef,truncated_idx,phiVal,A,b,case_num, T, N_xi_inUse)
% get function values from coefficeint of basis (evaluated at mesh)
K_est_val                = phiVal*coef;
structName.K_est_val     = K_est_val;
structName.c_hat         = coef;
structName.truncated_idx = truncated_idx;
if nargin>4
    f_hat     = reshape(A*coef, case_num, T, N_xi_inUse);
    temp_diff = (A*coef-b);
    structName.resnorm  = sum(temp_diff.^2);
    structName.residual = temp_diff;
    structName.f_hat    = f_hat;
end
end