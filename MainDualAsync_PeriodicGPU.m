%% define I/O
[os, ~, ~] = computer;
if strcmp(os,'MACI64')
    Projdir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/STDP_Project';
    Gitdir = '~/STDP';
elseif strcmp(os,'GLNXA64')
    % Projdir = '/gpfs/data/glimcherlab/BoShen/STDP_Project';
    % Gitdir = '/gpfs/data/glimcherlab/BoShen/STDP';
    Projdir = '/home/bs3667/STDP';
    Gitdir = '/home/bs3667/STDP';
elseif strcmp(os, 'PCWIN64')
    Projdir = 'C:\Users\Bo\Dropbox (NYU Langone Health)\Bo Shen Working files\STDP_Project';
    Gitdir = 'C:\Users\Bo\Documents\GitHub\STDP';
end
gnrloutdir = fullfile(Projdir, 'General');
Svmat_dir = fullfile(Gitdir, 'Simulations');
addpath(genpath(Gitdir));% Test main
ProjectName = 'DualASyncTestGPU';
%% Define time vector
dt = .1; % ms, time precision for simulation, in unit of second
duration = 51000; % ms
% Time vector
time = [0:dt:duration]';
timesteps = numel(time);
%% Spike train of the input, example trials
Ntrial = 6;
rng(29);
Seq = CreateEvents(Ntrial, dt/1000);
Seq = Seq(1:timesteps, 1:2)';

%% runner
Runner;

