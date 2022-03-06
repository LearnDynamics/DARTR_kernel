
function [Est_vector] = vector_kernel_run2(kernelInfo, regressionData, SAVE_DIR, plotOn,saveON,Bmat_type,normType,plot4paper)
% learn kernel by regression with regularization
if ~exist('saveON','var'); saveON=0; end


%% 1. learn kernel by regression 
[Est_vector, h_fig_coef] = vector_learn_kernel(regressionData, plotOn,Bmat_type,normType);  % learn kernel by regression 

%% 2 compute L2 error of estimator 
if isfield(kernelInfo, 'K_true')
    K_true   = kernelInfo.K_true; 
    rseq        = Est_vector.r_seq; % dx*(1:inferInfo.bdry_width)
    dr          = rseq(2)-rseq(1);  
    K_true_val  = K_true(rseq);
    if plotOn ==1 
        figure(h_fig_coef); hold on; 
        plot(rseq,K_true_val,'k:', 'Displayname', 'True', 'Linewidth', 3); 
    end
    Est_vector.K_true_val = K_true_val;

    fprintf('\n Vector Estimator error: \n ');   % plot4paper=1; 
    [~,~,L2_errors,h]= vector_getAll_L2_error2(Est_vector.cLcurveAll,Est_vector.rho_val,K_true_val',dr,plotOn,plot4paper);
    Est_vector.L2_errors = L2_errors; 
    
    if plot4paper ==1;  figure(h); figname = [SAVE_DIR,'figures/Est_vec_',kernelInfo.kernel_type]; 
        set_positionFontsAll_narrow1; legend 'boxoff'; % legend 'off'
        print([figname,'.pdf'],'-dpdf', '-bestfit');
    end
    
    if saveON==1
        filename = [SAVE_DIR,'vector_Est_',regressionData.data_str,'.mat'];
        save(filename,'Est_vector','regressionData','kernelInfo','L2_errors');
    end
end  

%% 3. compute loss values of each estimator
lossVal = getLoss_values_vec(Est_vector); 
Est_vector.lossval= lossVal; 
if plotOn ==1
    fprintf('\n Vector Estimator loss values: '); disp(lossVal);
end
end

