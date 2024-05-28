%% define I/O
DefineIO;
%% Define time vector
dt = 1; % ms, time precision for simulation, in unit of second
duration = 2068001; % ms
% Time vector
time = [dt:dt:duration]';
timesteps = numel(time);
ProjectName = sprintf('SavingTestGPU_%1.1fs', duration/1000);
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
%% Setup for visualization etc
Setup;
if gpuDeviceCount > 0
    gpuparallel = 1;
else
    gpuparallel = 0;
end
%% Build the neural network
show = 1;
if gpuparallel
    NetworkgeneratorPeriodicGPU;
else
    NetworkgeneratorPeriodic;
end

%% Specify project name and output
plotdir = fullfile(Projdir, ProjectName);
if ~exist(plotdir, 'dir')
    mkdir(plotdir);
end
Rsltfile = fullfile(plotdir,'Rslts.mat');

%% Inputs of single source, assume this long-range projection only intervenes Exct neurons
%% Define input connection matrix
Ntwk.Input.Source = 2; % Number of input source(s)
Ntwk.Input.N = 400*Ntwk.Input.Source; % number of the input (excitarory) neurons
Ntwk.Input.Tube = 100; % um, the diameter of the long-range projection range, shaped as a tube.
Ntwk.Input.AxonRange = 150; % um, the assumed standard deviation of axon physical connection range of a single input neuron
% assume the long-range projection target the center field of the tested patch
gpurng(2024);
r = gpuArray.randn(Ntwk.Input.N,1)*Ntwk.Input.Tube; % radius of the input fiber
gpurng(2025);
phi = gpuArray.rand(Ntwk.Input.N,1)*2*pi; % angle of the location of the input fiber
xE = cos(phi).*r;
shift = ([0:1/(1+Ntwk.Input.Source):1]-1/2)*Ntwk.Scale*2;
xE = xE + repmat(shift(2:end-1)', Ntwk.Input.N/Ntwk.Input.Source,1);
yE = sin(phi).*r;
Origins = repmat([1;2], Ntwk.Input.N/Ntwk.Input.Source, 1);
[xE, I] = sort(xE);
yE = yE(I);
Origins = Origins(I);
Ntwk.Input.Location = [xE, yE];
Ntwk.Input.Origins = Origins;
clear xE yE Origins r phi I;
% connections to the local excitatory neurons

[XInput, XE] = meshgrid(Ntwk.Input.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent local and columns represent Input projection
[YInput, YE] = meshgrid(Ntwk.Input.Location(:,2), Ntwk.Exct.Location(:,2));
DstcInput = sqrt(min(Ntwk.Scale*2 - abs(XInput - XE), abs(XInput - XE)).^2 + min(Ntwk.Scale - abs(YInput - YE), abs(YInput - YE)).^2); % Euclidean distance between each pair of neurons
p_Input = exp(-.5*(DstcInput/Ntwk.Input.AxonRange).^2); % probability of physical connection based on distance
gpurng(2024);
Ntwk.Cnnct_Input = p_Input >= gpuArray.rand(size(p_Input)); % projections from input, 0 or 1
clear XInput YInput XE YE DstcInput p_Input;
if show
    Ntwk = InputTuning(Ntwk, OKeeffe, gnrloutdir);
end
%% Visualize input sequences
h = figure;
filename = 'InputDynamic';
subplot(3,1,1); hold on;
for ii = Ntwk.Input.Source:-1:1
    plot(time'/1000, Seq(ii,:)*.9+(ii), '-', "Color", OKeeffe(ii,:), 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0 duration/1000 .5, Ntwk.Input.Source*1.45]); % Adjust the axis for better visualization
% ylim([.5, Ntwk.Input.Source*1.45]);
yticks([1:Ntwk.Input.Source]);
mysavefig(h, "Dummy Input", plotdir, 12, [3,6], 1);
% Simulate multiple neurons under driven of the same input signal
% Parameters for Poisson distribution
spikeRate = 100/1000; % spikes per milisecond
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = spikeRate * dt;

% Poisson generator: generate random numbers and compare to spike probability
leftt = floor(min(find(Seq(1,:)>0, 1 ), find(Seq(2,:)>0, 1 )*dt)/100)*100;
rightt = ceil((max(find(Seq(1,:)>0, 1 ), find(Seq(2,:)>0, 1 ))*dt+500)/100)*100;
tmpsteps = (rightt - leftt)/dt;
InputSpikes = gpuArray.rand(Ntwk.Input.N, tmpsteps);
trunck = 50000;
for i = 1:ceil(tmpsteps/trunck)
    timevec = 1+(i-1)*trunck:min(i*trunck, tmpsteps);
    InputSpikes(Ntwk.Input.Origins == 1, timevec) = InputSpikes(Ntwk.Input.Origins == 1, timevec) < spikeProbability*Seq(1, timevec + leftt/dt);
    InputSpikes(Ntwk.Input.Origins == 2, timevec) = InputSpikes(Ntwk.Input.Origins == 2, timevec) < spikeProbability*Seq(2, timevec + leftt/dt);
end

% Plot the raster plot
subplot(3,1,2);
hold on;
for n = 1:5:Ntwk.Input.N
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(n, :));
    % Convert indices to time
    spikeTimes = time(spikeIndices+leftt/dt);
    % Plot spikes for this neuron
    if Ntwk.Input.Origins(n) == 1
        plot(spikeTimes, n * ones(size(spikeTimes)), '.', "Color", OKeeffe(1,:), 'MarkerSize', 2);
    elseif Ntwk.Input.Origins(n) == 2
        plot(spikeTimes, n * ones(size(spikeTimes)), '.', "Color", OKeeffe(2,:), 'MarkerSize', 2);
    end
end
hold off;
xlabel('Time (ms)');
ylabel('Input neurons');
title('Spike train');
axis([leftt rightt 0.5 Ntwk.Input.N+0.5]); % Adjust the axis for better visualization
mysavefig(h, "Dummy Input", plotdir, 12, [3,6], 1);
clear InputSpikes spikeTimes spikeIndices SumStrength;
%% Simulating the neural network
% Initializing the network status at t=0
% - membrane potentials
ExctV = Ntwk.VL*gpuArray.rand(Ntwk.Exct.N,1);
InhbtV = Ntwk.VL*gpuArray.rand(Ntwk.Inhbt.N,1);
ExctRefraction = gpuArray.rand(Ntwk.Exct.N,1);
InhbtRefraction = gpuArray.rand(Ntwk.Inhbt.N,1);
% - synaptic conductance (presynaptic-activity dependent)
ExctgE = gpuArray.rand(Ntwk.Exct.N, 1);
ExctgI = gpuArray.rand(Ntwk.Exct.N, 1);
InhbtgE = gpuArray.rand(Ntwk.Inhbt.N, 1);
InhbtgI = gpuArray.rand(Ntwk.Inhbt.N, 1);
% - synaptic weights (plastic according to STDP rules)
WEI = Ntwk.wEI_initial; %rand(Ntwk.Inhbt.N, Ntwk.Exct.N);
WIE = Ntwk.wIE_initial; %rand(Ntwk.Exct.N, Ntwk.Inhbt.N);
WEE = Ntwk.wEE_initial; %rand(Ntwk.Exct.N, Ntwk.Exct.N);
% - spiking events
Espikes = gpuArray.rand(Ntwk.Exct.N,1);
Ispikes = gpuArray.rand(Ntwk.Inhbt.N,1);
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
% - intermediate variable for STPD convolution over time
xEpre = gpuArray.rand(1,Ntwk.Exct.N);
xEpost = gpuArray.rand(Ntwk.Exct.N,1);
xIpre = gpuArray.rand(1,Ntwk.Inhbt.N);
xIpost = gpuArray.rand(Ntwk.Inhbt.N,1);
% Dynamic variables of the example neurons
exampleE = Ntwk.Smpl.E; % E1, E2, EShare1, EShare2
exampleI = Ntwk.Smpl.I; % I1, I2, IShare
Exmpl.ExctV = [ExctV(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.InhbtV = [InhbtV(exampleI,1)'; gpuArray.rand(timesteps-1, numel(exampleI))];
Exmpl.ExctgE = [ExctgE(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.ExctgI = [ExctgI(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.InhbtgE = [InhbtgE(exampleI,1)'; gpuArray.rand(timesteps-1, numel(exampleI))];
Exmpl.InhbtgI = [InhbtgI(exampleI,1)'; gpuArray.rand(timesteps-1, numel(exampleI))];
Exmpl.WEI = gpuArray.rand(timesteps, numel(exampleI), numel(exampleE));
Exmpl.WEI(1,:,:) = WEI(exampleI,exampleE);
Exmpl.WIE = gpuArray.rand(timesteps, numel(exampleE), numel(exampleI));
Exmpl.WIE(1,:,:) = WIE(exampleE,exampleI);
Exmpl.WEE = gpuArray.rand(timesteps, numel(exampleE), numel(exampleE));
Exmpl.WEE(1,:,:) = WEE(exampleE,exampleE);
Exmpl.Espikes = [Espikes(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.Ispikes = [Ispikes(exampleI,1)'; gpuArray.rand(timesteps-1, numel(exampleI))];
Exmpl.xEpre = [xEpre(1,exampleE); gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.xEpost = [xEpost(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))];
Exmpl.xIpre = [xIpre(1,exampleI); gpuArray.rand(timesteps-1, numel(exampleI))];
Exmpl.xIpost = [xIpost(exampleI,1)'; gpuArray.rand(timesteps-1, numel(exampleI))];

% Other parameters
Ib = 150; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
ExctNoise = Ib + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
InhbtNosie = Ib + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
Exmpl.ExctNoise = [ExctNoise(exampleE,1)'; gpuArray.rand(timesteps-1, numel(exampleE))]; % ExctNoise(exampleE,1)
% simulation start...

ExctVec = 1:Ntwk.Exct.N;
InhbtVec = [1:Ntwk.Inhbt.N]*4;
for t = 1:(timesteps-1)
    % input spikes
%     InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
%     InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(1,t);
%     InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(2,t);
end
%% Save results
close all;
clearvars -except 'Ntwk' 'Seq' 'Exmpl' 'WEI' 'WIE' 'WEE' 'Rsltfile';
save(Rsltfile, 'Ntwk', 'Seq', 'Exmpl', 'WEI', 'WIE', 'WEE', '-v7.3');
