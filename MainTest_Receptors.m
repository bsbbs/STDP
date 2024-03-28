% Test main
%% define I/O
[os, ~, ~] = computer;
if strcmp(os,'MACI64')
    Projdir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/STDP_Project';
    Gitdir = '~/STDP';
elseif strcmp(os,'GLNXA64')
    Projdir = '/gpfs/data/glimcherlab/BoShen/Noise';
    Gitdir = '/gpfs/data/glimcherlab/BoShen/Noise';
end
gnrloutdir = Projdir;
Svmat_dir = fullfile(Gitdir, 'Simulations');
addpath(genpath(Gitdir));
%% setup
Setup;

%% Build network of single neurons
show = 1;
Networkgenerator;

%% define LIF model
Ntwk.VL = -70; % mV, resting potential
Ntwk.Vth = -50; % mV, threshold potential (Wang 2002; Vogels et al., 2011)
Ntwk.Vreset = -55; % mV, reset potential (Wang 2002)
Ntwk.VE = 0; % mV, excitatory synaptic potential
Ntwk.VI = -80; % mV, inhibitory synaptic potential (Vogels et al. 2011; Burkitt et al., 2004)
Ntwk.Exct.Cm = .5; % nF, membrane capacity for pyramidal neurons (Wang 2002)
Ntwk.Inhbt.Cm = .2; % nF, membrane capacity for interneurons (Wang 2002)
Ntwk.Exct.gL = 25; % nS, membrane leaky conductance for pyramidal neurons (Wang 2002)
Ntwk.Exct.gAMPA = .05; % nS, Wang 2002
Ntwk.Exct.gNMDA = .165; % nS, Wang 2002
Ntwk.Exct.gGABA = 1.3; % nS, Wang 2002
Ntwk.Inhbt.gL = 20; % nS, membrane leaky conductance for interneurons (Wang 2002)
Ntwk.Inhbt.gAMPA = .04; % nS, Wang 2002
Ntwk.Inhbt.gNMDA = .13; % nS, Wang 2002
Ntwk.Inhbt.gGABA = 1.0; % nS, Wang 2002
Ntwk.Exct.taum = Ntwk.Exct.Cm/Ntwk.Exct.gL*1000; % 20 ms, membrane time constant of pyramidal neurons (Wang, 2002; Vogels et al., 2011; Burkitt et al., 2004)
Ntwk.Inhbt.taum = Ntwk.Inhbt.Cm/Ntwk.Inhbt.gL*1000; % 10 ms, membrane time constant of interneurons (Wang 2002)
Ntwk.AMPA.tau = 2; % ms, the decay time constant of AMPA currents. (Wang 2002)
Ntwk.NMDA.tau.decay = 100; % ms, the decay time constant of NMDA currents. (Wang 2002)
Ntwk.NMDA.tau.rise = 2; % ms, the rising time constant of NMDA currents. (Wang 2002)
Ntwk.NMDA.alpha = .5; % ms^-1, rising speed of NMDA
Ntwk.GABA.tau = 5; % ms, the decay time constent of GABA currents. (Wang 2002)
Ntwk.Exct.tauREF = 2; % ms, refractory period for pyramidal neurons (Wang 2002)
Ntwk.Inhbt.tauREF = 1; % ms, refractory period for interneurons (Wang 2002)
%% define plasticity kernels
% Kernel, STDP sliding

%% Simulator
plotdir = fullfile(gnrloutdir, 'InitialTest');
if ~exist(plotdir, 'dir')
    mkdir(plotdir);
end
dt = 1; % ms, time precision for simulation, in unit of second
duration = 5000; % ms
% Time vector
time = [0:dt:duration]';
timesteps = numel(time);
%% Dummy inputs
numNeurons = 200;
% Parameters for Poisson distribution
spikeRate = 5/1000; % spikes per milisecond
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = spikeRate * dt;

% Generate random numbers and compare to spike probability
InputSpikes = rand(numNeurons,length(time)) < spikeProbability;

% Plot the raster plot
h = figure;
hold on;
for n = 1:numNeurons
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(n,:));
    % Convert indices to time
    spikeTimes = time(spikeIndices);
    % Plot spikes for this neuron
    plot(spikeTimes, n * ones(size(spikeTimes)), 'k.', 'MarkerSize', 5);
end
hold off;
xlabel('Time (ms)');
ylabel('Neuron');
title('Spike Trains of the Input');
axis([0 duration 0.5 numNeurons+0.5]); % Adjust the axis for better visualization
mysavefig(h, "Dummy Input", plotdir, 12, [3,2], 1);

%% single neuron fires
sNMDA = 0; % synaptic activity of NMDA
x = 0; % NMDA rising phase
sAMPA = 0; % synaptic activity of AMPA
sGABA = 0; % synaptic activity of GABA
w = 1;
refractionPeriod = Ntwk.Exct.tauREF/dt;
V = nan(timesteps,1);
V(1) = Ntwk.VL;
INMDA = zeros(timesteps,1);
IAMPA = zeros(timesteps,1);
IGABA = zeros(timesteps,1);
Isyn = zeros(timesteps,1);
spikes = zeros(timesteps,1);

last_spike_time = -Inf;
for t = 1:(timesteps-1)
    % NMDA receptor, synaptic conductance dynamics
    dx = -x/Ntwk.NMDA.tau.rise*dt + sum(InputSpikes(:,t));
    x = x + dx;
    dsNMDA = (-sNMDA/Ntwk.NMDA.tau.decay + Ntwk.NMDA.alpha*x*(1 - sNMDA))*dt;
    sNMDA = sNMDA + dsNMDA; 
    INMDA(t) = Ntwk.Exct.gNMDA*(Ntwk.VE - V(t))*sum(w*sNMDA);
    % AMPA receptor, synaptic conductance dynamics
    dsAMPA = -sAMPA/Ntwk.AMPA.tau*dt + sum(InputSpikes(:,t));
    sAMPA = sAMPA + dsAMPA; 
    IAMPA(t) = Ntwk.Exct.gAMPA*(Ntwk.VE - V(t))*sum(w*sAMPA);
    % dynamics of synaptic conductance of GABA receptor
    dsGABA = (-sGABA/Ntwk.GABA.tau + sum(InputSpikes(:,t)))*dt;
    sGABA = sGABA + dsGABA;
    IGABA(t) = Ntwk.Exct.gGABA*(Ntwk.VI - V(t))*sum(sGABA);
    % Overall current
    Isyn(t) = INMDA(t) + IAMPA(t) + IGABA(t);
    % Membrane potential change
    if (t - last_spike_time) > refractionPeriod
        dV = (Ntwk.Exct.gL*(Ntwk.VL - V(t)) + Isyn(t))*dt/Ntwk.Exct.Cm;
        V(t+1) = V(t) + dV;
    else
        V(t+1) = Ntwk.Vreset;
    end
    % fire
    if V(t+1) > Ntwk.Vth
        spikes(t+1) = 1;
        last_spike_time = t+1;
    end
  
end
%
h = figure;
subplot(2,1,1);
plot(time, V);
subplot(2,1,2);
plot(time, INMDA);