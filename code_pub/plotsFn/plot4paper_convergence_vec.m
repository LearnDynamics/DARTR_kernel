function plot4paper_convergence_vec(loss_L2err, noise_seq, dx_seq, kernel_type, reg_type, figname)
%{
loss_L2err: a cell with structures, each = 
    L2error_Ind: 2
          L2err: [2×3×5×5 double]
            res: [1×3×5×5 double]
%}
Nkernels = numel(kernel_type);
Res    = cell(Nkernels,1);
L2err  = cell(Nkernels,1);


f1 = @(st) squeeze(st.res);
f2 = @(st) squeeze(st.L2err(st.L2error_Ind,:,:,:));
Res   = cellfun(f1, loss_L2err, 'UniformOutput', false);
L2err = cellfun(f2, loss_L2err, 'UniformOutput', false);

[Nest, Ndx, Nnsr] = size(Res{1});

% plot settings
if Nest == 4;    space_subplots = 0.02;  
elseif Nest <=3; space_subplots = 0.02;  
end
ind_noise = [2, 4];    % select the index of noise to plot:  all if not plotting loss
cell_min    = @(st) min(st(:,:,ind_noise),[], 'all')*0.9;
cell_max    = @(st) max(st(:,:,ind_noise),[], 'all')*1.1;
ylim_lb_Res = cellfun(cell_min, Res, 'UniformOutput', false); 
ylim_ub_Res = cellfun(cell_max, Res, 'UniformOutput', false); 
ylim_lb_L2err = cellfun(cell_min, L2err, 'UniformOutput', false); 
ylim_ub_L2err = cellfun(cell_max, L2err, 'UniformOutput', false); 
n_noise_inplot = length(ind_noise); 
linestyle_color_left = {'-o','-.o'}; 
linestyle_color_right= {'-d','-.d'}; 

figure;
ha = tight_subplot(Nkernels,Nest,space_subplots,[0.1,0.08],[0.09,0.06]);   % figname = [fig_path,'funEst_all',name_str]; 
idx = 1;
lgndNoise = cell(n_noise_inplot,2); %1: L2err; 2: loss
for row = 1:Nkernels
    for col = 1:Nest
        axes(ha(idx)); hold on;  
        if row == 1, title(reg_type{col});end
        if  row == Nkernels %  && col == floor((Nest)/2)
            xlabel('\Delta x =0.0125\times\{1,2,4,8,16\}','Unit', 'normalized', 'FontSize', 18); % 'Position',[.2 0],
        end
        if col == 1 
            pos_subplot = get(ha(idx),'position'); lowerleftcorner = pos_subplot(2); h = pos_subplot(4);
            pos = lowerleftcorner+0.5*h+0.08*row^(1/4); 
            annotation('textarrow',[0.5 0.5],[0.5 0.5],'string', kernel_type{row}, ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[0.0001 pos 0 0],'FontSize',12,'FontWeight','bold');
        end
        for n = 1:n_noise_inplot
            n_noise = ind_noise(n);
            lgndNoise{n,1} = sprintf('nsr = %1.1f, error',noise_seq(n_noise));
            lgndNoise{n,2} = sprintf('nsr = %1.1f, loss',noise_seq(n_noise));
            yyaxis left; yticklabels('auto');    % yticks(10.^[-5:2]); 
            loglog(dx_seq, L2err{row}(col,:,n_noise),linestyle_color_left{n},'linewidth',1.5, 'Displayname', lgndNoise{n,1}); 
            xlim([dx_seq(1)*0.9, dx_seq(end)*1.1])
            ylim([ylim_lb_L2err{row} ylim_ub_L2err{row}]); grid on;
            yyaxis right;   % yticks(10.^[-5:2]); 
            loglog(dx_seq, Res{row}(col,:,n_noise),linestyle_color_right{n},'linewidth',1.5, 'Displayname', lgndNoise{n,2} );  
            ylim([ylim_lb_Res{row} ylim_ub_Res{row}]); 
        end
        idx = idx+1;
        yyaxis left;   set(gca, 'YScale', 'log');    set(gca, 'XScale', 'log');
        if col == 1; yticklabels('auto');  ylabel('L2 error'); else; yticklabels(''); end
        yyaxis right;  set(gca, 'YScale', 'log'); set(gca, 'XScale', 'log');
        if col == Nest;  yticklabels('auto'); ylabel('Loss value'); else; yticklabels(''); end
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








