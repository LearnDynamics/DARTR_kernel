 figname = 'lcurve_rkhs_mfOpt_Gaussian0125'; %%% l2, LL2, rkhs
 savefig(figname)


%% set position, size of figure; and set fonts
%  function set_positionFontsAll(nrow,ncol)
ftsz = 16; %mksz = 12; % if save to -dpdf
% ftsz =12; mksz = 8; % if save to -depsc
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
%set(findall(gcf,'-property','MarkerSize'),'MarkerSize',mksz);

%% get the number of subplots in gcf; 
cells_axes = findobj(gcf,'type','axes');
N = length(cells_axes);
pos1 = zeros(N,1); pos2 = zeros(N,1);
for n = 1:N
    pos1(n) = cells_axes(n).Position(1);
    pos2(n) = cells_axes(n).Position(2);
end
Ncols = numel(unique(pos1));
Nrows = numel(unique(pos2));

%% set figure size
width = 400+100*(Ncols-1); height = 400+200*(Nrows-1); 
set(gcf, 'Position',  [100, 1000, width, height]);  % set figure position + size [x , x, width height]
%% save figure
%pbaspect([1 1 1])
tightfig; 
if exist('figpath', 'var') && exist('figname', 'var')
    print([figpath,figname,'.pdf'],'-dpdf', '-bestfit'); 
end
if (~exist('figpath', 'var')) && exist('figname', 'var')
    print([figname,'.pdf'],'-dpdf', '-fillpage');  
end