% Test main
%% define I/O
[os, ~, ~] = computer;
if strcmp(os,'MACI64')
    Projdir = '/Users/bs3667/Dropbox (NYU Langone Health)/Bo Shen Working files/STDP_Project';
    Gitdir = '~/STDP';
elseif strcmp(os,'GLNXA64')
    Projdir = '/gpfs/data/glimcherlab/BoShen/STDP_Project';
    Gitdir = '/gpfs/data/glimcherlab/BoShen/STDP';
elseif strcmp(os, 'PCWIN64')
    Projdir = 'C:\Users\Bo\Dropbox (NYU Langone Health)\Bo Shen Working files\STDP_Project';
    Gitdir = 'C:\Users\Bo\Documents\GitHub\STDP';
end
gnrloutdir = Projdir;
Svmat_dir = fullfile(Gitdir, 'Simulations');
addpath(genpath(Gitdir));
%% Setup for visualization etc
Setup;

%% Build the neural network
show = 1;
Networkgenerator;

%% Specify project name and output
plotdir = fullfile(gnrloutdir, 'SingleInputTest');
if ~exist(plotdir, 'dir')
    mkdir(plotdir);
end

%% Define time vector
dt = .1; % ms, time precision for simulation, in unit of second
duration = 5000; % ms
% Time vector
time = [0:dt:duration]';
timesteps = numel(time);

%% Inputs of single source, assume this long-range projection only intervenes Exct neurons
%% Define input connection matrix
Ntwk.Input.Source = 1; % Number of input source(s)
Ntwk.Input.N = 800; % number of the input (excitarory) neurons
Ntwk.Input.Tube = 100; % um, the diameter of the long-range projection range, shaped as a tube.
Ntwk.Input.AxonRange = 150; % um, the assumed standard deviation of axon physical connection range of a single input neuron
% assume the long-range projection target the center field of the tested patch
rng(2024);
xE = randn(Ntwk.Input.N,1)*Ntwk.Input.Tube;
[xE, I] = sort(xE);
yE = rand(Ntwk.Input.N,1)*100 - 50; % zeros(Ntwk.Exct.N,1); % assume all neurons approximately on the same layer
Ntwk.Input.Location = [xE, yE];
clear xE yE;
% connections to the local excitatory neurons

[XInput, XE] = meshgrid(Ntwk.Input.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent local and columns represent Input projection
[YInput, YE] = meshgrid(Ntwk.Input.Location(:,2), Ntwk.Exct.Location(:,2));
DstcInput = sqrt((XInput - XE).^2 + (YInput - YE).^2); % Euclidean distance between each pair of neurons
p_Input = exp(-.5*(DstcInput/Ntwk.Input.AxonRange).^2); % probability of physical connection based on distance
Ntwk.Cnnct_Input = p_Input >= rand(size(p_Input)); % projections from input, 0 or 1 
clear XInput YInput XE YE DstcInput p_Input;
if show
    h = figure;
    filename = 'InputTuning';
    subplot(2,1,1); hold on;
    % connection strength distribution from each input
    legendLabels = cell(1, Ntwk.Input.Source);
    for i = 1:Ntwk.Input.Source
        SumStrength = sum(Ntwk.Cnnct_Input,2);
        plot(Ntwk.Exct.Location(:,1), SumStrength/max(SumStrength), 'LineWidth', 2, 'Color', OKeeffe(i,:));
        legendLabels{i} = sprintf('Input %d', i);
    end
    legend(legendLabels, 'Location', 'best');
    xlabel('\mum');
    ylabel('Coupling strength');
    mysavefig(h, filename, plotdir, 12, [3, 3], 2);
    
    subplot(2,1,2); hold on;
    undersample = 4;
    plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'kv', 'MarkerSize', 3);
    plot(Ntwk.Exct.Location(exampleE,1), Ntwk.Exct.Location(exampleE,2), 'cv', 'MarkerSize', 3);
    for i = 1:Ntwk.Input.Source
        for j = 1:undersample:Ntwk.Input.N 
            plot([Ntwk.Input.Location(j,1), Ntwk.Input.Location(j,1)], [100, Ntwk.Input.Location(j,2)], 'Color', OKeeffe(i,:), 'LineWidth',.5);
        end
    end
    xlabel('\mum');
    ylabel('\mum');
    xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
    ylim([-50, 100]);
    mysavefig(h, filename, plotdir, 12, [3, 3], 2);
end

%% Spike train of the input, example trials
Ntrial = 3;
rng(2024);
Seq = CreateEvents(Ntrial, dt/1000);
Seq = Seq(1:timesteps, 1)';
h = figure;
filename = 'InputDynamic';
subplot(3,1,1);
for ii = 1:Ntwk.Input.Source
    plot(time', Seq, '-', "Color", OKeeffe(ii,:), 'LineWidth',2);
end
title('Input signal');
xlabel('Time (ms)');
ylabel('Input');
ylim([-.5, Ntwk.Input.Source*1.5]);
mysavefig(h, "Dummy Input", plotdir, 12, [3,6], 1);
% Simulate multiple neurons under driven of the same input signal
% Parameters for Poisson distribution
spikeRate = 100/1000; % spikes per milisecond
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = spikeRate * dt;

% Poisson generator: generate random numbers and compare to spike probability
InputSpikes = rand(Ntwk.Input.N, timesteps) < spikeProbability*Seq;

% Plot the raster plot
subplot(3,1,2);
hold on;
for n = 1:50:Ntwk.Input.N
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(n, :));
    % Convert indices to time
    spikeTimes = time(spikeIndices);
    % Plot spikes for this neuron
    plot(spikeTimes, n * ones(size(spikeTimes)), 'k.', 'MarkerSize', 2);
end
hold off;
xlabel('Time (ms)');
ylabel('Neuron');
title('Spike train');
axis([0 duration 0.5 Ntwk.Input.N+0.5]); % Adjust the axis for better visualization
mysavefig(h, "Dummy Input", plotdir, 12, [3,6], 1);

%% Simulating the neural network
% Initializing the network status at t=0
% - membrane potentials
ExctV = Ntwk.VL*ones(Ntwk.Exct.N,1);
InhbtV = Ntwk.VL*ones(Ntwk.Inhbt.N,1);
ExctRefraction = zeros(Ntwk.Exct.N,1);
InhbtRefraction = zeros(Ntwk.Inhbt.N,1);
% - synaptic conductance (presynaptic-activity dependent)
ExctgE = zeros(Ntwk.Exct.N, 1);
ExctgI = zeros(Ntwk.Exct.N, 1);
InhbtgE = zeros(Ntwk.Inhbt.N, 1);
InhbtgI = zeros(Ntwk.Inhbt.N, 1);
% - synaptic weights (plastic according to STDP rules)
WEI = Ntwk.wEI_initial; %rand(Ntwk.Inhbt.N, Ntwk.Exct.N);
WIE = Ntwk.wIE_initial; %rand(Ntwk.Exct.N, Ntwk.Inhbt.N);
WEE = Ntwk.wEE_initial; %rand(Ntwk.Exct.N, Ntwk.Exct.N);
% - spiking events
Espikes = zeros(Ntwk.Exct.N,1);
Ispikes = zeros(Ntwk.Inhbt.N,1);
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
% - intermediate variable for STPD convolution over time
xEpre = zeros(1,Ntwk.Exct.N);
xEpost = zeros(Ntwk.Exct.N,1);
xIpre = zeros(1,Ntwk.Inhbt.N);
xIpost = zeros(Ntwk.Inhbt.N,1);
% Dynamic variables of the example neurons
Exmpl.ExctV = [ExctV(exampleE,1); zeros(timesteps-1, 1)];
Exmpl.InhbtV = [InhbtV(exampleI,1); zeros(timesteps-1, 1)];
Exmpl.ExctgE = [ExctgE(exampleE,1); zeros(timesteps-1, 1)];
Exmpl.ExctgI = [ExctgI(exampleE,1); zeros(timesteps-1, 1)];
Exmpl.InhbtgE = [InhbtgE(exampleI,1); zeros(timesteps-1, 1)];
Exmpl.InhbtgI = [InhbtgI(exampleI,1); zeros(timesteps-1, 1)];
Exmpl.WEI = [WEI(exampleI,exampleE); zeros(timesteps-1, 1)];
Exmpl.WIE = [WIE(exampleE,exampleI); zeros(timesteps-1, 1)];
Exmpl.WEE = [WEE(exampleE,exampleE); zeros(timesteps-1, 1)];
Exmpl.Espikes = [Espikes(exampleE,1); zeros(timesteps-1, 1)];
Exmpl.Ispikes = [Ispikes(exampleI,1); zeros(timesteps-1, 1)];
Exmpl.xEpre = [xEpre(1,exampleE); zeros(timesteps-1, 1)];
Exmpl.xEpost = [xEpost(exampleE,1); zeros(timesteps-1, 1)];
Exmpl.xIpre = [xIpre(1,exampleI); zeros(timesteps-1, 1)];
Exmpl.xIpost = [xIpost(exampleI,1); zeros(timesteps-1, 1)];

% Other parameters
Ib = 200; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
ExctNoise = Ib + randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
InhbtNosie = Ib + randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
Exmpl.ExctNoise = [ExctNoise(exampleE,1); zeros(timesteps-1, 1)]; % ExctNoise(exampleE,1)
% simulation start...
h = figure; hold on;
filename = 'RealtimeMonitor';
ExctVec = 1:Ntwk.Exct.N;
InhbtVec = [1:Ntwk.Inhbt.N]*4;
for t = 1:(timesteps-1)
    % Synaptic plasticity
    xEpre = xEpre -(xEpre/Ntwk.ExctSTDP.tau_prepost)*dt + Espikes';
    Exmpl.xEpre(t+1) = xEpre(exampleE);
    xEpost = xEpost -(xEpost/Ntwk.ExctSTDP.tau_postpre)*dt + Espikes;
    Exmpl.xEpost(t+1) = xEpost(exampleE);
    xIpre = xIpre -(xIpre/Ntwk.InhbtSTDP.tau_prepost)*dt + Ispikes';
    Exmpl.xIpre(t+1) = xIpre(exampleI);
    xIpost = xIpost -(xIpost/Ntwk.InhbtSTDP.tau_postpre)*dt + Ispikes;
    Exmpl.xIpost(t+1) = xIpost(exampleI);
    WEI = WEI + Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xIpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Ispikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Ispikes*xEpre + Ntwk.ExctSTDP.intercept_post*Ispikes*ones(size(Espikes'))); % pre -> post
    Exmpl.WEI(t+1) = WEI(exampleI, exampleE);
    WIE = WIE + Ntwk.InhbtSTDP.eta*(Ntwk.InhbtSTDP.sign_postpre*xEpost*Ispikes' + Ntwk.InhbtSTDP.intercept_pre*ones(size(Espikes))*Ispikes' ... % post -> pre
        + Ntwk.InhbtSTDP.sign_prepost*Espikes*xIpre + Ntwk.InhbtSTDP.intercept_post*Espikes*ones(size(Ispikes'))); % pre -> post
    Exmpl.WIE(t+1) = WIE(exampleE, exampleI);
    WEE = WEE + Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xEpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Espikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Espikes*xEpre + Ntwk.ExctSTDP.intercept_post*Espikes*ones(size(Espikes'))); % pre -> post
    Exmpl.WEE(t+1) = WEE(exampleE+1, exampleE);
    
    % Synaptic activities
    ExctgE = ExctgE - ExctgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Exct neurons
    + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_Input*InputSpikes(:,t) + Ntwk.Cnnct_EE.*WEE*Espikes);
    Exmpl.ExctgE(t+1) = ExctgE(exampleE);
    ExctgI = ExctgI - ExctgI/Ntwk.Synapse.tauInhbt*dt ... % inhibitory synaptic conductance on Exct neurons
    + Ntwk.Synapse.gbarI*(Ntwk.Cnnct_IE.*WIE*Ispikes);
    Exmpl.ExctgI(t+1) = ExctgI(exampleE);
    InhbtgE = InhbtgE - InhbtgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Inhbt neurons
    + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_EI.*WEI*Espikes);
    Exmpl.InhbtgE(t+1) = InhbtgE(exampleI);
    InhbtgI = InhbtgI - InhbtgI/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Inhbt neurons
    Exmpl.InhbtgI(t+1) = InhbtgI(exampleI);
    
    % Updating OU noise
    ExctNoise = ExctNoise + ((Ib - ExctNoise)/Ntwk.Noise.tauN + randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm)*dt;
    Exmpl.ExctNoise(t+1) = ExctNoise(exampleE);
    InhbtNosie = InhbtNosie + ((Ib - InhbtNosie)/Ntwk.Noise.tauN + randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm)*dt;
    % Membrane potential change for Exct neurons
    ExctRefraction = ExctRefraction - 1;
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + ExctgE.*(Ntwk.VE - ExctV) + ExctgI.*(Ntwk.VI - ExctV) + ExctNoise)/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    ExctV(Espikes) = Ntwk.Vfire;
    Exmpl.ExctV(t+1) = ExctV(exampleE);
    Exmpl.Espikes(t+1) = Espikes(exampleE);
    ExctRefraction(Espikes) = refractionPeriod.E;
    plot(time(t)*ones(sum(Espikes),1), ExctVec(Espikes), 'k.', 'MarkerSize', 2);

    % Membrane potential change for Inhbt neurons
    InhbtRefraction = InhbtRefraction - 1;
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + InhbtgE.*(Ntwk.VE - InhbtV) + InhbtgI.*(Ntwk.VI - InhbtV) + InhbtNosie)/Ntwk.Inhbt.Cm*dt;
    InhbtV = InhbtV + dV;
    InhbtV(InhbtRefraction>0) = Ntwk.Vreset;
    
    % Inhbt neurons fire
    Ispikes = InhbtV > Ntwk.Vth;
    InhbtV(Ispikes) = Ntwk.Vfire-10;
    Exmpl.InhbtV(t+1) = InhbtV(exampleI);
    Exmpl.Ispikes(t+1) = Ispikes(exampleI);
    InhbtRefraction(Ispikes) = refractionPeriod.I;
    plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 2);
end
axis([1700 2310 0 Ntwk.Exct.N]);
xlabel('Time (ms)');
ylabel('Neurons');
title('Raster plot of Exct/Inhbt neurons');
mysavefig(h, filename, plotdir, 12, [8, 6], 1);

%% visualizing example neurons
h = figure;
filename = 'Example neurons activity';
subplot(2,1,1); hold on;
lg = [];
lg(1) = plot(time, Exmpl.ExctV, 'k-');
% plot(time(Exmpl.Espikes==1), Ntwk.Vfire*Exmpl.Espikes(Exmpl.Espikes==1), 'k.', 'MarkerSize',4);
lg(2) = plot(time, Exmpl.InhbtV, 'r-');
% plot(time(Exmpl.Ispikes==1), Ntwk.Vfire*Exmpl.Ispikes(Exmpl.Ispikes==1), 'r.', 'MarkerSize',4);
legend(lg, {'Exct', 'Inhbt'}, 'Location','eastoutside');
% axis([1700 2400 -80 11]);
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Activity of a single neuron');
mysavefig(h, filename, plotdir, 12, [8, 4], 1);

subplot(2,1,2); hold on;
lg = [];
lg(1) = plot(time, Exmpl.ExctgE, 'k-');
lg(2) = plot(time, Exmpl.ExctgI, 'k--');
lg(3) = plot(time, Exmpl.InhbtgE, 'r-');
lg(4) = plot(time, Exmpl.InhbtgI, 'r--');
legend(lg, {'gE on Exct cell','gI on Exct cell', 'gE on Inhbt cell','gI on Inhbt cell'}, 'Location','eastoutside');
% axis([1700 2400 0 20]);
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Conductance (nS)');
title('Synaptic activity');
mysavefig(h, filename, plotdir, 12, [8,4], 1);

h = figure;
filename = 'Example synaptic weights';
subplot(4,1,1); hold on;
lg = [];
lg(1) = plot(time, Exmpl.xEpre, 'k-');
lg(2) = plot(time, Exmpl.xEpost, 'k--');
lg(3) = plot(time, Exmpl.xIpre, 'r-');
lg(4) = plot(time, Exmpl.xIpost, 'r--');
legend(lg, {'xEpre','xEpost', 'xIpre','xIpost'}, 'Location','eastoutside');
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Kernel of weight');
title('STDP Integration');
mysavefig(h, filename, plotdir, 12, [5,8], 1);
subplot(4,1,2); hold on;
plot(time, Exmpl.WEI, 'k-');
legend('WEI', 'Location','eastoutside');
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Plastic weight');
mysavefig(h, filename, plotdir, 12, [5,8], 1);
subplot(4,1,3); hold on;
plot(time, Exmpl.WIE, 'r-');
legend('WIE', 'Location','eastoutside');
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Plastic weight');
mysavefig(h, filename, plotdir, 12, [5,8], 1);
subplot(4,1,4); hold on;
plot(time, Exmpl.WEE, 'k--');
legend('WEE', 'Location','eastoutside');
xlim([1700 2400]);
xlabel('Time (ms)');
ylabel('Plastic weight');
mysavefig(h, filename, plotdir, 12, [5,8], 1);

%% visualize weight change
h = figure;
filename = 'WEI_change_5s';
imagesc(WEI - Ntwk.wEI_initial)
colormap(bluewhitered);
colorbar;
xlabel("Exct neurons");
ylabel("Inhbt neurons");
mysavefig(h, filename, plotdir, 12, [5,4], 1);

h = figure;
filename = 'WIE_change_5s';
imagesc(WIE - Ntwk.wIE_initial)
colormap(bluewhitered);
colorbar;
xlabel("Inhbt neurons");
ylabel("Exct neurons");
mysavefig(h, filename, plotdir, 12, [5,4], 1);

h = figure;
filename = 'WEE_change_5s';
imagesc(WEE - Ntwk.wEE_initial)
colormap(bluewhitered);
colorbar;
xlabel("Exct neurons");
ylabel("Exct neurons");
mysavefig(h, filename, plotdir, 12, [5,4], 1);

