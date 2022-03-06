function [L2_err,L2rho_err,L2_errors,h]= vector_getAll_L2_error2(cLcurveAll,rho_val,Ktrue,dx,plotOn,plot4paper)
%  get function value from coefficients and basis values 
%%

lebesgueWt   = 1/(dx*length(rho_val)); 
linestyle_color = {'--','-.','--','-.','--',':',':',':',':',':'}; 
h = [];
% cLcurveAll = Est.cLcurveAll;
% field_str_all = cLcurveAll{'creg_rkhs1','creg_rkhs2','creg_rkhs3','creg_H1'};
name_str     = cLcurveAll.normType;   % {'rkhs1','rkhs2','rkhs3','H1'};

Nstr         = length(name_str); 
L2_err       = zeros(1,Nstr);  L2rho_err = 0*L2_err; 
Kest_array   = zeros(length(rho_val),Nstr); 

if plotOn ==1
    xi = dx* (1:length(rho_val));
    h = figure; plot(xi,Ktrue,'k:','linewidth',2);  hold on;
end
for nn = 1:Nstr
    field_str     = ['creg_',name_str{nn}];
    Kest          = cLcurveAll.(field_str); % K_est_val_truncated_l2;
    L2_err(nn)    = sqrt(sum((Ktrue - Kest).^2 * dx));
    L2rho_err(nn) = sqrt(sum((Ktrue - Kest).^2 .*rho_val * dx));
    Kest_array(:,nn) = Kest;
    if plotOn==1
        figure(h); plot(xi,Kest,linestyle_color{nn},'linewidth',2);
    end
end
if plotOn ==1
    figure(h);  legend('True',name_str{1:end},'Location','best');
    xlabel('r'); % ylabel('Kernel value');
    diff = max(Ktrue)-min(Ktrue);
    axis([0,xi(end),min(Ktrue)-diff*0.2,max(Ktrue)+diff*0.2])
    if plot4paper ==1
        normalizer = max(Ktrue)/2/max(rho_val); rho_length = length(rho_val);
        nx = length(xi); gap = ceil(nx/rho_length); indx = 1:gap:nx; xirho = xi(indx); % sparser index
        % area(xirho,rho_val(indx)*normalizer, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeAlpha', 0,'Displayname', '\rho');
        area(xirho,rho_val(indx)*normalizer, 'FaceColor', 'c', 'FaceAlpha', 0.2,'EdgeAlpha', 0,'Displayname', '\rho');
        yyaxis right; yticklabels('auto');
    end
end
L2_errors.L2_err    = L2_err*sqrt(lebesgueWt);
L2_errors.L2rho_err = L2rho_err;
L2_errors.name_str  = name_str;
disp(name_str);
disp([L2_err; L2rho_err]);

end

 
    


