function resArray = getLoss_values_vec(Est)
% compute the loss values of estimators for vector estimator. For fucntion
% estimator, need to select optimal dimension. 
%{

%}

% Est_dx_seq{1}.Est_vector{1,1}.L2_errors.name_str
%    {'Lcurve-l2'}    {'Lcurve-L2'}    {'Lcurve-RKHS'}    {'svd-l2'}    {'svd-L2'}    {'svd-RKHS'}
name_str     = Est.L2_errors.name_str;
N_est        = numel(name_str);
resArray     = zeros(1,N_est);    % entry are ordered the same as the above name_str; 
 
Abar          = Est.matInfo.Abar;
bbar          = Est.matInfo.bbar;
b_norm_sq     = Est.matInfo.b_norm_sq;

if isfield(Est, 'cLcurveAll')
    for i = 1:N_est
        cLcurveAll = Est.cLcurveAll;
        s = split(name_str{i}, '_');
        field_name = ['creg_', s{end}];
        resArray(i) = computeLoss(cLcurveAll.(field_name),Abar,bbar,b_norm_sq);
    end
end

end