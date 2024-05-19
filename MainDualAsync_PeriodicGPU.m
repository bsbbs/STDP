% Test main
ProjectName = 'DualASyncTestGPU';
%% Define time vector
dt = .1; % ms, time precision for simulation, in unit of second
duration = 51000; % ms
% Time vector
time = [0:dt:duration]';
timesteps = numel(time);
%% Spike train of the input, example trials
Ntrial = 6;
rng(29);
Seq = CreateEvents(Ntrial, dt/1000);
Seq = Seq(1:timesteps, 1:2)';

%% runner
Runner;

