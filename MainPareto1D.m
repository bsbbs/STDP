%% define I/O
DefineIO1D;

%% Build the neural network
visualize = 1;
NetworkgeneratorGPU1D;
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
Testdir = fullfile(Projdir, 'ParetoValues');
if ~exist(Testdir, 'dir')
    mkdir(Testdir);
end

%% Parrellel simulations
Covs = 0:.1:1;
evs = [evspool(:,1) evspool(:,1)];
Type = 'Pareto';
Inputstruct = Ntwk.Input;
mypool = parpool(11);
parfor runi = 1:11
    % define values
    Cov = Covs(runi);
    Runningdir = fullfile(Testdir, sprintf('%s_Cov%1.1f', Type, Cov));
    if ~exist(Runningdir,'dir')
        mkdir(Runningdir);
    end
    fprintf('%s, Val = %1.1f\n', Type, Cov);
    valfile = fullfile(Runningdir, sprintf('Seq_%s_Cov%1.1f_%itrials.mat', Type, Cov, Ntrial));
    values = ParetoSequence(Ntrial, Cov);
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
    RunnerNMDA1D(Ntwk, Seq, dt, Runningdir);
    % RunnerNMDA1D;
    fprintf('%s, Cov = %1.1f Done\n', Type, Cov);
    s = struct("WEE", WEE, "WEI", WEI, "WIE", WIE);
    save(Trainedfile, '-fromstruct', s);
end
