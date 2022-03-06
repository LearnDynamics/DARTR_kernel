
figure;
%subplot(131);
plot(Est.xi, Est.RKHSbasis_vals(:,1:5),'linewidth', 1.5); title('1-5')
legend('\psi_1','\psi_2','\psi_3','\psi_4','\psi_5', 'Location', 'best','Orientation','horizontal', 'NumColumns', 2)
%subplot(132);
%plot(Est.xi, Est.RKHSbasis_vals(:,26:27),'linewidth', 1.5); title('26-27')
%legend('\psi_{26}','\psi_{27}', 'Location', 'best','Orientation','horizontal')  
%subplot(133);
%plot(Est.xi, Est.RKHSbasis_vals(:,51:52),'linewidth', 1.5); title('51-52')
%legend('\psi_{51}','\psi_{52}', 'Location', 'best','Orientation','horizontal');
figpath = [SAVE_DIR,'figures/'];
figname = ['RKHSbasisAdaptive', '_ui_', obsInfo.u_str, '_xj_', obsInfo.x_mesh_str, '_H_', inferInfo.basis_str];
set_positionFontsAll;


%{
%% psi at x=0
figure;
semilogy(abs(Est.RKHSbasis_vals(1,:)), 'o-', 'Linewidth', 1.5); hold on;
semilogy(abs(Est_final.RKHSbasis_vals(1,:)), 'd-','Linewidth', 1.5); hold off;
ylabel(' |\psi_i(0.05)|'); xlabel('i');
title('RKHS basis');
legend('Uniform', 'Adaptive');
figpath = [SAVE_DIR,'figures/'];
figname = ['RKHSbasisAdaptive', '_ui_', obsInfo.u_str, '_xj_', obsInfo.x_mesh_str, '_H_', inferInfo.basis_str];
set_positionFontsAll;
%}