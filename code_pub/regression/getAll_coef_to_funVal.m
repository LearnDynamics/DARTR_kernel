function [EstVal]= getAll_coef_to_funVal(Est,phiVal)
%  get function value from coefficients and basis values
%{
c_plain        = Abar\bbar;           % ***** good when no noise
c_regAll.c_plain        = c_plain;
c_regAll.cSVD_l2_max     = cSVD_l2_max;      % trucated SVD: ratio < threshold
c_regAll.cSVD_L2_max     = cSVD_L2_max;      % ***** for ill-posed inverse: regularization when b has noise/error > smallest eigenvalue
c_regAll.cSVD_rkhs_max   = cSVD_rkhs_max;
c_regAll.cSVD_all        = cSVDratio_all; % a structure with L2-RKHS estimator information
c_regAll.cLcurve_rkhs  = cLcurveAll.creg_rkhs;         % ***** for ill-posed inverse: regularization when b has noise/error > smallest eigenvalue
c_regAll.cLcurve_l2    = cLcurveAll.creg_l2;
c_regAll.cLcurve_L2    = cLcurveAll.creg_L2;
%}

%% Just for 3RKHS and H1
if isfield(Est, 'cLcurveAll')
    normType_list = {'l2', 'L2', 'rkhs1','rkhs2','rkhs3','RKHS', 'H1', 'H2', 'gRKHS'};
    cLcurveAll =Est.cLcurveAll;
    for i = 1:numel(normType_list)
        s = normType_list{i};
        creg_field = ['creg_', s];
        K_est_val_field = ['K_est_val_Lcurve_', s];
        if isfield(cLcurveAll, creg_field)
        coef  = cLcurveAll.(creg_field);
        [Kest,~] = coef_to_funVal(coef,0,phiVal);
        EstVal.(K_est_val_field) = Kest;
        end
    end
    
end
%% for whole reg pack
if isfield(Est, 'c_regAll')
    c_regAll = Est.c_regAll;
    
    % truncated SVD
    truncated_idx = [];
    coef  = c_regAll.cSVD_l2_max;
    [Kest_truncated_l2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cSVD_L2_max;
    [Kest_truncated_L2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cSVD_H1_max;
    [Kest_truncated_H1,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cSVD_H2_max;
    [Kest_truncated_H2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cSVD_rkhs_max;
    [Kest_truncated_rkhs,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    
    % Lcurve
    coef  = c_regAll.cLcurve_l2;
    [Kest_Lcurve_l2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cLcurve_L2;
    [Kest_Lcurve_L2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cLcurve_rkhs;
    [Kest_Lcurve_rkhs,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cLcurve_H1;
    [Kest_Lcurve_H1,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    coef  = c_regAll.cLcurve_H2;
    [Kest_Lcurve_H2,~] = coef_to_funVal(coef,truncated_idx,phiVal);
    
    
    
    EstVal.K_est_val_Lcurve_l2      = Kest_Lcurve_l2;
    EstVal.K_est_val_Lcurve_L2      = Kest_Lcurve_L2;
    EstVal.K_est_val_Lcurve_RKHS    = Kest_Lcurve_rkhs;
    EstVal.K_est_val_Lcurve_H1    = Kest_Lcurve_H1;
    EstVal.K_est_val_Lcurve_H2    = Kest_Lcurve_H2;
    
    EstVal.K_est_val_truncated_l2   = Kest_truncated_l2;
    EstVal.K_est_val_truncated_L2   = Kest_truncated_L2;
    EstVal.K_est_val_truncated_H1   = Kest_truncated_H1;
    EstVal.K_est_val_truncated_H2   = Kest_truncated_H2;
    EstVal.K_est_val_truncated_RKHS = Kest_truncated_rkhs;
end
end