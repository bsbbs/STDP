%% Get the SLURM array task ID
taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
fprintf('Running MATLAB TASK %d\n', taskID);

%% define I/O
DefineIO1D;

%% Build the neural network
visualize = 1;
Ntwk = NetworkgeneratorGPU1D(gnrloutdir, visualize, Mycolors);
close all;
%% Define input sequences
Ntrial = 1600; % numbers of testing trials
evfile = fullfile(gnrloutdir, sprintf('Events%itrials.mat', Ntrial));
if ~exist(evfile, 'file')
    evs = Generator(Ntrial);
    save(evfile, 'evs');
else
    load(evfile, 'evs');
end

%% Define values of the input events
Testdir = fullfile(Projdir, 'ConstantValues');
if ~exist(Testdir, 'dir')
    mkdir(Testdir);
end

%% Parrellel simulations
ValAmps = repmat([.01, .05, .1, .5, 1, 2, 8, 30, 50, 100],1, 2);
evsSync = [evs(:,1) evs(:,1)];
evsAsync = evs;
runi = taskID;
sessi = ceil(runi/10);
if sessi == 1
    Type = 'Sync';
    Seq = evsSync;
elseif sessi == 2
    Type = 'Async';
    Seq = evsAsync;
end
% define values
ValAmp = ValAmps(runi);
Runningdir = fullfile(Testdir, sprintf('%s_ValAmp%3.3f', Type, ValAmp));
if ~exist(Runningdir,'dir')
    mkdir(Runningdir);
end
fprintf('%s, ValAmp = %3.3f\n', Type, ValAmp);
valfile = fullfile(Runningdir, sprintf('%s_ValAmp%3.3f_%itrials.mat', Type, ValAmp, Ntrial));
if ~exist(valfile, 'file')
    % sigma = 0;
    % values = ParetoSequence(Ntrial, sigma);
    values = ones(Ntrial,2)*ValAmp;
    save(valfile, 'values');
end
% Visualize input sequences, integrated with values
%     if visualize
%         VisualizeSeq(Seq, values, Mycolors, Runningdir);
%     end

%     %% runner
%     dt = .1; % ms, time precision for simulation
%     duration = length(SeqPool)*dt; % ms
%     time = (dt:dt:duration)';
%     timesteps = numel(time);
%     save(Seqfile, 'SeqPool','evs','values','time','duration','timesteps');
%     Ntwk = RunnerNMDA1D(Ntwk, Seq, dt, Runningdir, visualize);
%     %% Evaluation
%     ChecktheTestNMDA1D;
%     ComputationNMDA1D;


