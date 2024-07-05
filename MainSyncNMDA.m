%% define I/O
DefineIO;
%% Define time vector, check input sequences
dt = .1; % ms, time precision for simulation, in unit of second
% Spike train of the input, example trials
Ntrial = 1600;
ProjectName = sprintf('SyncNMDA_%i', Ntrial);
plotdir = fullfile(Projdir, ProjectName);
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Seqfile = fullfile(plotdir, 'Seq.mat');
if ~exist(Seqfile, 'file')
    % sigma = 0;
    % values = ParetoSequence(Ntrial, sigma);
    values = ones(Ntrial,2);
    [Seq, evs] = Generator(Ntrial, values, dt);
    Seq = Seq(:, 1:2);
    % Time vector
    duration = length(Seq)*dt; % ms
    time = [dt:dt:duration]';
    timesteps = numel(time);
    save(Seqfile, 'Seq','evs','values','time','duration','timesteps');
else
    load(Seqfile);
end
duration = length(Seq)*dt; % ms
h = figure;
filename = 'InputDynamic';
subplot(3,1,1); hold on;
for ii = 2:-1:1
    plot(time/1000, Seq(:,ii)*.9+(ii), '-', 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0, time(end)/1000, .5, 2*1.45]); % Adjust the axis for better visualization
% ylim([.5, Ntwk.Input.Source*1.45]);
yticks([1:2]);
%close(h);
%% runner
RunnerNMDA;

