
function plot4paper_l2L2rkhs(figname,Nnoise,dx_seq,L2errEst_ftn,linestyle_color,RegType,noise_ratios,lossval)
N_est = numel(L2errEst_ftn(1,:,1,1));
figure;
space_subplots = 0.05;  
ind_noise = 1:Nnoise;    % select the index of noise to plot:  all if not plotting loss
if exist('lossval','var') 
    space_subplots = 0.10; 
    ind_noise = [1,Nnoise-1];   % select the index of noise to plot: here only the first ant the 2nd-to-last noise when plotting loss
    ylim_lb_loss = min(lossval(1,1:N_est,:,ind_noise),[], 'all')*0.9; 
    ylim_ub_loss = max(lossval(1,1:N_est,:,ind_noise),[],'all')*1.1; 
end
ha = tight_subplot(1,N_est,space_subplots);   % figname = [fig_path,'funEst_all',name_str]; 
n_noise_inplot = length(ind_noise); 
lgndNoise = cell(1,n_noise_inplot);  %   

ylim_lb = min(L2errEst_ftn(1,1:N_est,:,ind_noise),[], 'all')*0.9; 
ylim_ub = max(L2errEst_ftn(1,1:N_est,:,ind_noise),[], 'all')*1.1; 
for ind = 1:N_est     % ind= 1,2,3 for Lcurve, 4-6 for SVD 
    axes(ha(ind));      title(RegType{ind}); hold on; 
    for n_noise =1:n_noise_inplot
        ind_inplot = ind_noise(n_noise);
        % plot the L2 error
        if exist('lossval','var');  yyaxis left;  end 
        err_temp = squeeze(L2errEst_ftn(1,ind,:,ind_inplot));
        loglog(dx_seq, err_temp,linestyle_color{ind_inplot},'linewidth',1.5 );  hold on;
        ax = gca;
        xticklabels('auto');  yticklabels('auto'); ax.XScale = 'log';ax.YScale = 'log';
        ylim([ylim_lb,ylim_ub]);   % manually set the ylim for vector estimators
        if ind ==1; lgndNoise{n_noise} = sprintf('nsr %2.1f',noise_ratios(ind_inplot)); end
%         if ind ==N_est
%             x     = log(dx_seq); y = log(err_temp);
%             [a,~] = polyfit(x(1:end),y(1:end),1);
%             lgndNoise{n_noise} = sprintf('nsr %2.1f, slope=%1.2f',noise_ratios(ind_inplot),a(1));
%         end
        
        % plot the loss function
        if exist('lossval','var');  yyaxis right          
            lossval_temp = squeeze(lossval(1,ind,:,ind_inplot));
            loglog(dx_seq, lossval_temp,linestyle_color{ind_inplot},'linewidth',1.5 );  hold on;
            ax = gca; yticklabels('auto'); ax.YScale = 'log';
            ylim([ylim_lb_loss,ylim_ub_loss]);   % manually set the ylim for vector estimators
        end
    end
    if ind ==1 %|| ind==N_est; 
        legend(lgndNoise,'location','northwest');  end 
    if exist('lossval','var') 
       if ind == 1; yyaxis left; ylabel('L2 error'); end
       if ind == N_est;  yyaxis right; ylabel('Loss value'); end
    else 
       if ind == 1 || ind==N_est; ylabel('L2 error'); end
    end
     xlabel('dx'); 
end
% set(ha(1:3),'XTickLabelMode','auto'); 
% set(ha(1:3),'YTickLabelMode','auto');% ha(2).YTickLabelMode = 'auto'; 
set_positionFontsAll_narrow; print([figname,'.pdf'],'-dpdf', '-bestfit'); % savefig(figname);

% set_positionFontsAll; % print([figname,'.eps'],'-depsc');
end

