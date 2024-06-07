%% define I/O
DefineIO;
%% Define time vector
dt = 1; % ms, time precision for simulation, in unit of second
duration = 4500000; % ms
% Time vector
time = [dt:dt:duration]';
timesteps = numel(time);
sigma = 0.1;
ProjectName = sprintf('DualPareto%1.1f_%1.0fs', sigma, duration/1000);
plotdir = fullfile(Projdir, ProjectName);
if ~exist(plotdir, 'dir')
    mkdir(plotdir);
end
% Spike train of the input, example trials
Ntrial = 1600;
[Seq, values] = ParetoSequences(Ntrial, sigma, dt/1000);
Seq = Seq(1:timesteps, 1:2)';
Seq(:,1) = Seq(:,1)/mean(values(:,1));
Seq(:,2) = Seq(:,2)/mean(values(:,2));
h = figure;
filename = 'InputDynamic';
subplot(4,1,1); hold on;
for ii = 2:-1:1
    plot(time'/1000, Seq(ii,:)*.09+(ii), '-', 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0 duration/1000 .5, 2*1.45]); % Adjust the axis for better visualization
yticks([1:2]);
mysavefig(h, filename, plotdir, 12, [3,6], 1);
subplot(4,2,3); hold on;
histogram(values(:,1)/mean(values(:,1)), 100);
subplot(4,2,4); hold on;
histogram(values(:,2)/mean(values(:,1)), 100);
subplot(2,1,2); hold on;
plot(Seq(1,:), Seq(2,:), '.');
xlabel('Channel 1');
ylabel('Channel 2');
mysavefig(h, filename, plotdir, 12, [3,6], 1);
%% runner
Runner;

