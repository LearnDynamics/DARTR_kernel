function [minRes, minRes_ind, Res,dimH_seq] = min_res_Est_ftn(Est)
NdimH = length(Est); 
name_str = Est{1}.L2_errors.name_str;
N_est = numel(name_str);
dimH_seq = zeros(NdimH, 1);
Res = zeros(N_est,NdimH); 

%name_ref;  

for m=1:NdimH
    for nn = 1:N_est
        if isfield(Est{m}, 'cLcurveAll')
            cLcurveAll = Est{m}.cLcurveAll;  
            s = split(name_str{nn}, '_'); 
            field_name = ['creg_', s{end}];
            c_est = cLcurveAll.(field_name);
        end 
        Abar = Est{m}.matInfo.Abar;
        bbar = Est{m}.matInfo.bbar;
        b_norm_sq = Est{m}.matInfo.b_norm_sq;
        Res(nn,m) = computeLoss(c_est,Abar,bbar,b_norm_sq);
        dimH_seq(m) = numel(c_est);
    end
end


% regularizer_name_str_all = cell(1,N_est);
% for m=1:NdimH
%     for nn = 1:ftn_Nstr 
%         regularizer_name_str_all{nn} = ftn_field_str_coef_all{ftn_select_ind(nn)}{2}; 
%         c_est = getfield(temp{m}, ftn_field_str_coef_all{ftn_select_ind(nn)}{:}); 
%         Abar = temp{m}.matInfo.Abar;
%         bbar = temp{m}.matInfo.bbar;
%         b_norm_sq = temp{m}.matInfo.b_norm_sq;
%         Res(nn,m) = computeLoss(c_est,Abar,bbar,b_norm_sq);
%         
%         dimH_seq(m) = numel(c_est);
%     end
% end
[minRes,minRes_ind] = min(Res,[],2); 


end