
function h_fig_kernel_est =  plot_kernel_est(EstVal,Est,plotOn)
 % plot the estimated kernels in EstVal  
 h_fig_kernel_est =[]; 
if plotOn
    h_fig_kernel_est = figure;
    %subplot(111);
    xi = Est.r_seq;
    plot(xi, EstVal.K_est_val, 'c','Displayname', 'LSE', 'Linewidth', 1);hold on;
    plot(xi, EstVal.K_est_val_Lcurve_l2, '--','Displayname', 'L-curve l^2 norm', 'Linewidth', 2);
    plot(xi, EstVal.K_est_val_Lcurve_L2, '--','Displayname', 'L-curve L^2 norm', 'Linewidth', 2);
    plot(xi, EstVal.K_est_val_Lcurve_RKHS,'--', 'Displayname', 'L-curve RKHS norm', 'Linewidth', 2);
    plot(xi, EstVal.K_est_val_truncated_l2,'-.', 'Displayname', 'Truncated l^2', 'Linewidth', 2);
    plot(xi, EstVal.K_est_val_truncated_L2,'-.', 'Displayname', 'Truncated L^2', 'Linewidth', 2);
    plot(xi, EstVal.K_est_val_truncated_RKHS,'-.', 'Displayname', 'Truncated RKHS', 'Linewidth', 2);
    % % super-impose rho
    % normalizer = max(K_est_val)/2/max(rho.rho_val); rho_length = length(rho.rho_val);
    % nx = length(xi); gap = ceil(nx/rho_length); indx = 1:gap:nx; xirho = xi(indx); % sparser index
    % area(xirho,rho.rho_val(indx)*normalizer, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
    hold off;
    legend(); ylim([-50, 50])
    title('Estimators')
    

end
end 
