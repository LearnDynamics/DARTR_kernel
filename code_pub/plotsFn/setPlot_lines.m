% setPlots for three lines plots

% matlab RGB  darl colors
matcolor = zeros(6,3);
matcolor(1,:) = [0 0.4470 0.7410];        % 'darkblue'
matcolor(2,:) = [0.8500 0.3250 0.0980];   % 'darkorange'
matcolor(3,:) = [0.9290 0.6940 0.1250];   % 'darkyellow'
matcolor(4,:) = [0.4940 0.1840 0.5560];   % 'darkpurple'
matcolor(5,:) = [0.4660 0.6740 0.1880];   % 'darkgreen'
matcolor(6,:) = [0.3010 0.7450 0.9330];   % 'lightblue'
code_str = {'#0072BD', '#D95319','#EDB120',	'#7E2F8E', '#77AC30','#4DBEEE'};
name_str = {'darkblue','darkorange','darkyellow','darkpurple','darkgreen','lightblue'};
dark_color.matcolor; 
dark_color.code_str = code_str; 
dark_color.name_str = name_str;




linestyle  = {'-.';'--';':'}; 
colorArray = {'b';'r';'k'}; 
lineW      = {1;1;2}; 

h_axes = gca; h_plot = get(h_axes,'Children');
set(h_plot,{'LineStyle'},linestyle, {'Color'}, colorArray,{'LineWidth'},lineW);


color_tripe =[[0 0.4470 0.7410];]
linestyle_color = {'r:','b--'}; 
xi = 1:1:10;
figure; plot(xi,sin(xi),linestyle_color{1}); 


% colorArray = {'#D95319','#EDB120','#0072BD','#77AC30','#4DBEEE','cyan'}; 
% lineW      = [1,2,2,2,2,2,2]; 
     h_axes = gca; h_plot = get(h_axes,'Children');
     set(h_plot, {'Color'}, colorArray,{'LineWidth'},lineW);
     