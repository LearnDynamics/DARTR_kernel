function plot4paper_convergence(loss_L2err_info, noise_seq, dx_seq, kernel_type, reg_type, figname, selectby)
%{
loss_L2err_info{all kernels, certain example} = 
  struct with fields:
     L2error_Ind: 1: L2 Lebsgue error; 2: L2(rho) error
        minL2err: [2×Nest×Ndx×Nnsr double]
    res_minL2err: [1×Nest×Ndx×Nnsr double]
          minRes: [1×Nest×Ndx×Nnsr double]--->plot this by default
    L2err_minRes: [2×Nest×Ndx×Nnsr double]--->plot this by default
%}
Nkernels = numel(kernel_type);
Res = cell(Nkernels,1);
L2err = cell(Nkernels,1);

if ~exist('selectby', 'var') || strcmp(selectby, 'MinRes')
    f1 = @(st) squeeze(st.minRes);
    f2 = @(st) squeeze(st.L2err_minRes(st.L2error_Ind,:,:,:));
    Res   = cellfun(f1, loss_L2err_info, 'UniformOutput', false);
    L2err = cellfun(f2, loss_L2err_info, 'UniformOutput', false);
elseif strcmp(selectby, 'MinL2err')
    f1 = @(st) squeeze(st.res_minL2err);
    f2 = @(st) squeeze(st.minL2err(st.L2error_Ind,:,:,:));
    Res   = cellfun(f1, loss_L2err_info, 'UniformOutput', false);
    L2err = cellfun(f2, loss_L2err_info, 'UniformOutput', false);
end
[Nest, Ndx, Nnsr] = size(Res{1});

% plot settings
space_subplots = 0.05;  
ind_noise = [1, Nnsr-1];    % select the index of noise to plot:  all if not plotting loss
cell_min = @(st) min(st(:,:,ind_noise),[], 'all')*0.9;
cell_max = @(st) max(st(:,:,ind_noise),[], 'all')*1.1;
ylim_lb_Res = cellfun(cell_min, Res, 'UniformOutput', false); 
ylim_ub_Res = cellfun(cell_max, Res, 'UniformOutput', false); 
ylim_lb_L2err = cellfun(cell_min, L2err, 'UniformOutput', false); 
ylim_ub_L2err = cellfun(cell_max, L2err, 'UniformOutput', false); 


 n_noise_inplot = length(ind_noise); 
linestyle_color_left = {'-o','-.o'}; 
linestyle_color_right= {'-d','-.d'}; 
figure;
ha = tight_subplot(Nkernels,Nest,space_subplots);   % figname = [fig_path,'funEst_all',name_str]; 
idx = 1;
lgndNoise = cell(n_noise_inplot,2); %1: L2err; 2: loss
for row = 1:Nkernels
    for col = 1:Nest
        axes(ha(idx)); hold on;  
        if row == 1, title(reg_type{col});end
        if col == floor((Nest+1)/2) && row == Nkernels
            xlabel('\Delta x = 0.0125\times\{1, 2, 4, 8, 16\}','Unit', 'normalized','Position',[1.2 0], 'FontSize', 20);
        end
        if col == 1 
            pos_subplot = get(ha(idx),'position'); lowerleftcorner = pos_subplot(2); h = pos_subplot(4);
            pos = lowerleftcorner+0.5*h+0.08*row^(1/4); 
            annotation('textarrow',[0.5 0.5],[0.5 0.5],'string', kernel_type{row}, ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[0.0001 pos 0 0],'FontSize',12,'FontWeight','bold');
        end
        for n = 1:n_noise_inplot
            n_noise = ind_noise(n);
            lgndNoise{n,1} = sprintf('nsr = %2.0f, error',noise_seq(n_noise));
            lgndNoise{n,2} = sprintf('nsr = %2.0f, loss',noise_seq(n_noise));
            yyaxis left; yticklabels('auto');    % yticks(10.^[-5:2]); 
            loglog(dx_seq, L2err{row}(col,:,n_noise),linestyle_color_left{n},'linewidth',1.5, 'Displayname', lgndNoise{n,1});  
            ylim([ylim_lb_L2err{row} ylim_ub_L2err{row}]); grid on;
            yyaxis right;   % yticks(10.^[-5:2]); 
            loglog(dx_seq, Res{row}(col,:,n_noise),linestyle_color_right{n},'linewidth',1.5, 'Displayname', lgndNoise{n,2} );  
            ylim([ylim_lb_Res{row} ylim_ub_Res{row}]);
        end
        idx = idx+1;
        yyaxis left;   set(gca, 'YScale', 'log');    set(gca, 'XScale', 'log');
        if col == 1; yticklabels('auto');  ylabel('L2 error');end
        yyaxis right;  set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
        if col == Nest;  yticklabels('auto'); ylabel('Loss value');end
        hold off;
        if col == floor((Nest+1)/2) && row == Nkernels
            legend('Position',[0.53 0.01 0 0],'Orientation','vertical','NumColumns',2*n_noise_inplot, 'FontSize', 20)
            legend box off;
        end
    end   
end


set_positionFontsAll_narrow;
print([figname,'.pdf'],'-dpdf', '-bestfit');
end










