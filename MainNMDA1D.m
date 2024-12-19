%% define I/O
DefineIO1D;
% Setup for visualization etc
ColorPalette;

%% Define time vector, check input sequences
dt = .1; % ms, time precision for simulation, in unit of second
% dt = 1;
% Spike train of the input, example trials
Ntrial = 1600;
ValAmps = [.01, .05, .1, .5, 1, 2, 8, 30, 50, 100];
mypool = parpool(numel(ValAmps)*2);
parfor runi = 1:20
    sessi = ceil(runi/numel(ValAmps));
    if sessi == 1
        ProjectName = sprintf('Sync');
    elseif sessi == 2
        ProjectName = sprintf('Async');
    end
    plotdir = fullfile(Projdir, ProjectName);
    if ~exist(plotdir,'dir')
        mkdir(plotdir);
    end
    Seqfile = fullfile(plotdir, 'Seq.mat');
    %if ~exist(Seqfile, 'file')
        % sigma = 0;
        % values = ParetoSequence(Ntrial, sigma);
        values = ones(Ntrial,2);
        [Seq, evs] = Generator(Ntrial, values, dt);
        Seq = Seq(:, [1, 2]);
        % Time vector
        duration = length(Seq)*dt; % ms
        time = [dt:dt:duration]';
        timesteps = numel(time);
        %save(Seqfile, 'Seq','evs','values','time','duration','timesteps');
   %else
        %load(Seqfile);
    %end
    % duration = length(Seq)*dt; % ms
    ValAmp = ValAmps(runi);
    subplotdir = fullfile(plotdir, 'ValAmp%3.3',ValAmp);
    if ~exist(subplotdir,'dir')
        mkdir(subplotdir);
    end
    fprintf('%s, ValAmp = %3.3f\n', ProjectName, ValAmp);
    %% runner
    % [Ntwk, WEI, WIE, WEE] = RunnerNMDA1D(ValAmp, gnrloutdir, plotdir, subplotdir);
    %% Evaluation
    %ChecktheTestNMDA1D;
    %ComputationNMDA1D;
end

