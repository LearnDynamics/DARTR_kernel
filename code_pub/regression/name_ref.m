% vec
vec_select_ind   = [5 6 7  2 3 4];
vec_name_str_all     = {'LSE','svd-l2','svd-L2','svd-RKHS','Lcurve-l2','Lcurve-L2','Lcurve-RKHS'}; 
vec_field_str_all    = {'cSVD_l2_all.c_plain','cSVD_l2_max','cSVD_L2_max','cSVD_rkhs_max',...
                    'cLcurve_l2','cLcurve_L2','cLcurve_rkhs'}; 
vec_name_str     = vec_name_str_all(vec_select_ind); 
vec_Nstr    = length(vec_name_str);


% ftn
ftn_select_ind   = [4 5 6 1 2 3]; % select from 1-6; make Lcurve-RKHS being the 4th
ftn_name_str_all     = {'svd-l2','svd-L2','svd-RKHS','Lcurve-l2','Lcurve-L2','Lcurve-RKHS'}; 
ftn_field_str_all    = {'K_est_val_truncated_l2','K_est_val_truncated_L2','K_est_val_truncated_RKHS',...
                'K_est_val_Lcurve_l2','K_est_val_Lcurve_L2','K_est_val_Lcurve_RKHS'}; 
ftn_field_str_coef_all    = {'c_regAll.cSVD_l2_max','c_regAll.cSVD_L2_max','c_regAll.cSVD_rkhs_max',...
                'c_regAll.cLcurve_l2','c_regAll.cLcurve_L2','c_regAll.cLcurve_rkhs'}; 
ftn_field_str_coef_all = cellfun(@(x)regexp(x,'\.','split'), ftn_field_str_coef_all, 'UniformOutput',false);
ftn_name_str     = ftn_name_str_all(ftn_select_ind); 
ftn_Nstr         = length(ftn_name_str); 

%% for compare 3RKHS and H1
ftn_name_str_all     = {'rkhs1','rkhs2','rkhs3','H1'}; 
ftn_field_str_all    = {'K_est_val_Lcurve_rkhs1','K_est_val_Lcurve_rkhs2','K_est_val_Lcurve_rkhs3','K_est_val_Lcurve_H1'}; 
ftn_field_str_coef_all    = {'cLcurveAll.creg_rkhs1','cLcurveAll.creg_rkhs2','cLcurveAll.creg_rkhs3','cLcurveAll.creg_H1'}; 
ftn_field_str_coef_all = cellfun(@(x)regexp(x,'\.','split'), ftn_field_str_coef_all, 'UniformOutput',false);
ftn_name_str     = ftn_name_str_all; 
ftn_Nstr         = length(ftn_name_str); 

%{
Est.c_regAll.
            c_plain: [199×1 double]
         cSVD_l2_max: [199×1 double]
         cSVD_L2_max: [199×1 double]
       cSVD_rkhs_max: [199×1 double]
         cSVD_l2_all: [1×1 struct]
     cSVD_L2rkhs_all: [1×1 struct]
      cSVD_l2_cumsum: [199×1 double]
      cSVD_L2_cumsum: [199×1 double]
    cSVD_rkhs_cumsum: [199×1 double]
        cLcurve_rkhs: [199×1 double]
          cLcurve_l2: [199×1 double]
          cLcurve_L2: [199×1 double]
           lambda_l2: 1.9650e-07
           lambda_L2: 2.0342e-07
         lambda_rkhs: 2.5855e-04
      inds_lmabda_l2: [198×1 double]
      inds_lmabda_L2: [199×1 double]
    inds_lmabda_rkhs: [199×1 double]
%}