function h = plot4paper_Picard(eigA,eigL,ratio_rkhs,ratioA,ratioL)
% plot the Picard ratio 
    %str_picard  = {'Eig-vals L',  'Eig-vals A', 'Ratio l2','Ratio L2','Ratio RKHS'}; 
    str_picard  = {'Eigenvalues L',  'Eigenvalues A','Picard ratio l2','Picard ratio L2','Picard ratio RKHS'}; 

    h= figure;
    hold on;
    p1 = semilogy(abs(eigL), '-o','linewidth', 1, 'Markersize',1.2, 'Color', '#0072BD');
    set(p1, 'MarkerFaceColor', get(p1,'Color')); 
    p2=semilogy(abs(eigA), '-o','linewidth', 1, 'Markersize',1.2, 'Color', '#77AC30');
    set(p2, 'MarkerFaceColor', get(p2,'Color')); 
    p4=semilogy(abs(ratioA), '-d','linewidth', 1, 'Markersize',1.2, 'Color','#D95319');
    set(p4, 'MarkerFaceColor', get(p4,'Color')); 
    p3=semilogy(abs(ratioL), '-d','linewidth', 1, 'Markersize',1.2, 'Color', '#EDB120');
    set(p3, 'MarkerFaceColor', get(p3,'Color')); 
    p5=semilogy(abs(ratio_rkhs), '-d','linewidth', 1, 'Markersize',1.2, 'Color', '#7E2F8E');
    set(p5, 'MarkerFaceColor', get(p5,'Color')); 
    h2= yline(1, 'k:', 'linewidth', 1,'HandleVisibility','off');
   
    
    xlim([1 inf]);    
    %ylim([min(min(abs(eigA)), min(eigL))*0.9, max([abs(ratioL);abs(ratioA);abs(ratio_rkhs);1])*1.2]); 
    ylim([10^-10, 10^3]); 
    %yticklabels({})
    legend(str_picard,'Location', 'southwest', 'NumColumns', 2, 'Orientation', 'horizontal');legend boxoff;
    set(gca, 'YScale', 'log')
    hold off; 
end