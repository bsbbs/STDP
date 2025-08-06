[os, ~, ~] = computer;
if strcmp(os,'MACA64')
    Projdir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/STDP_Project/Analyses';
    Gitdir = '~/STDP';
elseif strcmp(os,'GLNXA64')
    % Projdir = '/gpfs/data/glimcherlab/BoShen/STDP_Project/Analyses';
    % Gitdir = '/gpfs/data/glimcherlab/BoShen/STDP';
    Projdir = '/scratch/bs3667/STDP_Project/Analyses';
    Gitdir = '/home/bs3667/STDP';
elseif strcmp(os, 'PCWIN64')
    Projdir = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project\Analyses';
    Gitdir = 'C:\Users\Bo\Documents\GitHub\STDP';
end
gnrloutdir = fullfile(Projdir, 'General1D');
if ~exist(gnrloutdir,'dir')
    mkdir(gnrloutdir);
end
addpath(genpath(Gitdir));
% check GPU
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    warning("Implementing this code requires GPU");
end
% loading color palette
ColorPalette;