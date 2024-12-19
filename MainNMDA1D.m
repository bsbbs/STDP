%% define I/O
DefineIO1D;
% Setup for visualization etc
ColorPalette;
%% Build the neural network
show = 1;
NetworkgeneratorPeriodicGPU1D;

%% Define time vector, check input sequences
dt = .1; % ms, time precision for simulation, in unit of second
% Spike train of the input, example trials
Ntrial = 1600;
Testdir = fullfile(Projdir, 'FixedValue');
if ~exist(Testdir, 'dir')
    mkdir(Testdir);
end
Seqfile = fullfile(Testdir, 'SeqPool.mat');
if ~exist(Seqfile, 'file')
    % sigma = 0;
    % values = ParetoSequence(Ntrial, sigma);
    values = ones(Ntrial,2);
    [SeqPool, evs] = Generator(Ntrial, values, dt);
    duration = length(SeqPool)*dt; % ms
    time = (dt:dt:duration)';
    timesteps = numel(time);
    save(Seqfile, 'SeqPool','evs','values','time','duration','timesteps');
else
    load(Seqfile);
end

%% Parrellel simulations
ValAmps = repmat([.01, .05, .1, .5, 1, 2, 8, 30, 50, 100],1, 2);
Seq1 = SeqPool(:, [1, 2]);
Seq2 = SeqPool(:, [1, 3]);
mypool = parpool(20);
parfor runi = 1:20
    sessi = ceil(runi/10);
    if sessi == 1
        ProjectName = 'Sync';
        Seq = Seq1;
    elseif sessi == 2
        ProjectName = 'Async';
        Seq = Seq2;
    end
    plotdir = fullfile(Testdir, ProjectName);
    if ~exist(plotdir,'dir')
        mkdir(plotdir);
    end
   
    ValAmp = ValAmps(runi);
    subplotdir = fullfile(plotdir, sprintf('ValAmp%3.3f', ValAmp));
    if ~exist(subplotdir,'dir')
        mkdir(subplotdir);
    end
    fprintf('%s, ValAmp = %3.3f\n', ProjectName, ValAmp);
    %% runner
    RunnerNMDA1D;
    %% Evaluation
    ChecktheTestNMDA1D;
    ComputationNMDA1D;
end

