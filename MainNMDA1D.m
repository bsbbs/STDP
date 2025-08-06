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
Testdir = fullfile(Projdir, 'ConstantValues');
if ~exist(Testdir, 'dir')
    mkdir(Testdir);
end

%% Parrellel simulations
Vals = repmat([.1, .5, 1, 2, 5, 10, 50], 1, 2);
evsSync = [evspool(:,1) evspool(:,1)];
evsAsync = evspool;
Inputstruct = Ntwk.Input;
mypool = parpool(14);
parfor runi = 1:14
    sessi = ceil(runi/7);
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
    RunnerNMDA1D(Ntwk, Seq, dt, Runningdir);
    % RunnerNMDA1D;
    fprintf('%s, Val = %3.2f Done\n', Type, Val);
    s = struct("WEE", WEE, "WEI", WEI, "WIE", WIE);
    save(Trainedfile, '-fromstruct', s);
end


%% Probe the compuational function of the network
% Vals = [.1, .5, 1, 2, 5, 10, 50];
Vals = [.1, 1, 5, 50];
Types = {'Sync','Async'};
for si = 1%:2
    for vi = 1:numel(Vals)
        Runningdir = fullfile(Testdir, sprintf('%s_Val%2.1f', Types{si}, Vals(vi)));
        fprintf('%s, Val = %2.1f\n', Types{si}, Vals(vi));
        fl = dir(fullfile(Runningdir, 'Weights_Seg*.mat'));
        numeric_parts_str = regexp({fl.name}, '\d+', 'match', 'once');
        file_numbers = str2double(numeric_parts_str);
        [~, max_index] = max(file_numbers);
        testf = fl(max_index);
        load(fullfile(testf.folder,testf.name));
        sim_dir = fullfile(Runningdir, sprintf('Test_Seg%i', max(file_numbers)));
        mkdir(sim_dir);
        ComputationNMDA1D;
    end
end

%% Explore possible mechanisms
Runningdir = fullfile(Testdir, 'Sync_Val1.0');
fl = dir(fullfile(Runningdir, 'Weights_Seg*.mat'));
numeric_parts_str = regexp({fl.name}, '\d+', 'match', 'once');
file_numbers = str2double(numeric_parts_str);
[~, max_index] = max(file_numbers);
testf = fl(max_index);
load(fullfile(testf.folder,testf.name));
sim_dir = fullfile(Runningdir, sprintf('Mechanism_Seg%i', max(file_numbers)));
mkdir(sim_dir);
