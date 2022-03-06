%% 
% 
% parentpath = pwd;
% 
% addpath([parentpath '/sysSolver/']);
% addpath([parentpath '/settings/']); 
% addpath([parentpath '/generateData/']);
% addpath([parentpath '/plotsFn/']);
% addpath([parentpath '/regularization/']);
% addpath([parentpath '/regression/']);
folder = fileparts(which(mfilename));
addpath(genpath(folder));

% addpath(genpath(fileparts(which(mfilename))));
clear folder
%% save data to your local DIR:  
% making DIR to your root DIR
if ispc
    SAVE_DIR = [getenv('USERPROFILE'), '\SIDA_RKHS\output'];
else
    SAVE_DIR = [getenv('HOME'),'/SIDA_RKHS/output/'];     
end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end    
if ~exist([SAVE_DIR,'figures/'],'dir'), mkdir([SAVE_DIR,'figures/']); end 
addpath(SAVE_DIR);


