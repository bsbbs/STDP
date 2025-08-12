%% define I/O
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
gnrloutdir = fullfile(Projdir, 'General_Isolated');
if ~exist(gnrloutdir,'dir')
    mkdir(gnrloutdir);
end
addpath(genpath(Gitdir));

% loading color palette
ColorPalette;%(gnrloutdir);
%% Build the neural network
visualize = 1;
Networkgenerator_Isolated;
close all;

%% Define input sequences
Ntrial = 1600; % numbers of testing trials
dt = .1; % ms, time precision for simulation
evfile = fullfile(gnrloutdir, sprintf('EventsPool%itrials.mat', Ntrial));
if ~exist(evfile, 'file')
    evspool = Generator(Ntrial);
    save(evfile, 'evspool');
else
    load(evfile, 'evspool');
end

%% Define values of the input events
Testdir = fullfile(Projdir, 'Isolated_ConstantValues');
if ~exist(Testdir, 'dir')
    mkdir(Testdir);
end

%% Parrellel simulations
Vals = repmat(1, 1, 2); % [.1, .5, 1, 2, 5, 10, 50]
evsSync = [evspool(:,1) evspool(:,1)];
evsAsync = evspool;
Inputstruct = Ntwk.Input;
% mypool = parpool(2);
for runi = 2
    sessi = runi; %ceil(runi/7);
    if sessi == 1
        Type = 'Sync';
        evs = evsSync;
    elseif sessi == 2
        Type = 'Async';
        evs = evsAsync;
    end
    % define values
    Val = Vals(runi);
    Runningdir = fullfile(Testdir, sprintf('%s_Val%2.1f', Type, Val));
    if ~exist(Runningdir,'dir')
        mkdir(Runningdir);
    end
    fprintf('%s, Val = %2.1f\n', Type, Val);
    valfile = fullfile(Runningdir, sprintf('Seq_%s_Val%2.1f_%itrials.mat', Type, Val, Ntrial));
    % sigma = 0;
    % values = ParetoSequence(Ntrial, sigma);
    values = ones(Ntrial,2)*Val;
    Seq = struct('evs',evs,'values',values);
    s = struct("Seq", Seq);
    save(valfile, '-fromstruct', s);

    % Visualize input sequences, integrated with values
    if visualize
        VisualizeSeq(Inputstruct, Seq, dt, Mycolors, Runningdir);
    end

    %% runner
    %% Specify project name and output
    Trainedfile = fullfile(Runningdir,'TrainedWeights.mat');
    [WEE, WEI, WIE] = RunnerIsolated(Ntwk, Seq, dt, Runningdir);
    % RunnerNMDA1D;
    fprintf('%s, Val = %3.2f Done\n', Type, Val);
    s = struct("WEE", WEE, "WEI", WEI, "WIE", WIE);
    save(Trainedfile, '-fromstruct', s);
end
