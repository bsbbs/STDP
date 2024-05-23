%% define I/O
DefineIO;
%% Define time vector
dt = 1; % ms, time precision for simulation, in unit of second
duration = 2068001; % ms
% Time vector
time = [dt:dt:duration]';
timesteps = numel(time);
ProjectName = sprintf('DualSyncTestGPU_%1.1fs', duration/1000);
% Spike train of the input, example trials
Ntrial = 258;
rng(29);
Seq = CreateEvents(Ntrial, dt/1000);
Seq = Seq(1:timesteps, 1:2)';
Seq(2,:) = Seq(1,:);
h = figure;
filename = 'InputDynamic';
subplot(3,1,1); hold on;
for ii = 2:-1:1
    plot(time'/1000, Seq(ii,:)*.9+(ii), '-', 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0 duration/1000 .5, 2*1.45]); % Adjust the axis for better visualization
% ylim([.5, Ntwk.Input.Source*1.45]);
yticks([1:2]);
%close(h);
%% runner
Runner;

