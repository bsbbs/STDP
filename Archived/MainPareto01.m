%% define I/O
DefineIO;
% Setup for visualization etc
ColorPalette;
%% Define time vector
dt = .1; % ms, time precision for simulation, in unit of second
% Spike train of the input, example trials
Ntrial = 1600;
sigma = 0.1;
ProjectName = sprintf('Pareto%1.1f_%i', sigma, Ntrial);
plotdir = fullfile(Projdir, ProjectName);
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Seqfile = fullfile(plotdir, 'Seq.mat');
if ~exist(Seqfile, 'file')
    rng(29);
    values = ParetoSequence(Ntrial, sigma);
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
h = figure;
filename = 'InputDynamic';
subplot(5,1,1); hold on;
for ii = 2:-1:1
    plot(time/1000, Seq(:,ii)*.2+(ii), '-', 'LineWidth', 1);
end
% title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([evs(110,1)-.5, evs(120,1)+1.5, .5, 2*1.45]); % Adjust the axis for better visualization
yticks([1:2]);
mysavefig(h, filename, plotdir, 12, [3,8], 1);
subplot(5,1,2); hold on;
for ii = 2:-1:1
    plot(time/1000, Seq(:,ii)*.09+(ii), '-', 'LineWidth', 1);
end
% title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0 duration/1000 .5, 2*1.45]); % Adjust the axis for better visualization
yticks([1:2]);
mysavefig(h, filename, plotdir, 12, [3,8], 1);
subplot(5,2,5); hold on;
histogram(values(:,1), 100);
xlabel('Value 1');
ylabel('Frequency');
mysavefig(h, filename, plotdir, 12, [3,8], 1);
subplot(5,2,6); hold on;
histogram(values(:,2), 100);
xlabel('Value 2');
ylabel('Frequency');
mysavefig(h, filename, plotdir, 12, [3,8], 1);
subplot(2,1,2); hold on;
scatter(values(:,1), values(:,2), 'filled');
xlabel('Channel 1');
ylabel('Channel 2');
mysavefig(h, filename, plotdir, 12, [3,8], 1);
%% runner
Runner;

