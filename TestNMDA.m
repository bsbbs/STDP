% Test NMDA synaptic conductance
%% initialize parameters
plotdir = fullfile(gnrloutdir, 'InitialTest');
dt = .1;
duration = 1000;
timevector = dt:dt:duration;
Nsteps = numel(timevector);
Ntwk.NMDA.tau.decay = 100; % ms, the decay time constant of NMDA currents. (Wang 2002)
Ntwk.NMDA.tau.rise = 2; % ms, the rising time constant of NMDA currents. (Wang 2002)
Ntwk.NMDA.alpha = .5; % ms^-1, rising speed of NMDA
%% define input signal
numNeurons = 200;
% Parameters for Poisson distribution
spikeRate = 5/1000; % spikes per milisecond
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = spikeRate * dt;

% Generate random numbers and compare to spike probability
InputSpikes = rand(numNeurons,length(timevector)) < spikeProbability;

InputSpikes(:,round(Nsteps/2):end) = 0;

h = figure;
subplot(3,1,1); hold on;
for n = 1:numNeurons
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(n,:));
    % Convert indices to time
    spikeTimes = timevector(spikeIndices);
    % Plot spikes for this neuron
    plot(spikeTimes, n * ones(size(spikeTimes)), 'k.', 'MarkerSize', 5);
end
hold off;
xlabel('Time (ms)');
ylabel('Neuron');
title('Spike Trains of the Input');
axis([0 duration 0.5 numNeurons+0.5]); % Adjust the axis for better visualization
mysavefig(h, "Test NMDA", plotdir, 12, [3,6], 1);
%% define synapse
x = zeros(Nsteps,1);
sNMDA = zeros(Nsteps,1);
for t = 1:(numel(timevector)-1)
    dx = -x(t)/Ntwk.NMDA.tau.rise*dt + sum(InputSpikes(:,t));
    x(t+1) = x(t) + dx;
    dsNMDA = (-sNMDA(t)/Ntwk.NMDA.tau.decay + Ntwk.NMDA.alpha*x(t)*(1 - sNMDA(t)))*dt;
    sNMDA(t+1) = sNMDA(t) + dsNMDA;
end
%% Visualize
subplot(3,1,2);
plot(timevector, x);
xlabel('Time (ms)');
ylabel('x');
title('Rising variable');
mysavefig(h, "Test NMDA", plotdir, 12, [3,6], 1);
subplot(3,1,3);
plot(timevector, sNMDA);
xlabel('Time (ms)');
ylabel('sNMDA');
title('Overall conductance');
mysavefig(h, "Test NMDA", plotdir, 12, [3,6], 1);