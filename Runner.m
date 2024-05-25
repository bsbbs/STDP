% Runner

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
ExctV = Ntwk.VL*gpuArray.ones(Ntwk.Exct.N,1);
InhbtV = Ntwk.VL*gpuArray.ones(Ntwk.Inhbt.N,1);
ExctRefraction = gpuArray.zeros(Ntwk.Exct.N,1);
InhbtRefraction = gpuArray.zeros(Ntwk.Inhbt.N,1);
% - synaptic conductance (presynaptic-activity dependent)
ExctgE = gpuArray.zeros(Ntwk.Exct.N, 1);
ExctgI = gpuArray.zeros(Ntwk.Exct.N, 1);
InhbtgE = gpuArray.zeros(Ntwk.Inhbt.N, 1);
InhbtgI = gpuArray.zeros(Ntwk.Inhbt.N, 1);
% - synaptic weights (plastic according to STDP rules)
WEI = Ntwk.wEI_initial; %rand(Ntwk.Inhbt.N, Ntwk.Exct.N);
WIE = Ntwk.wIE_initial; %rand(Ntwk.Exct.N, Ntwk.Inhbt.N);
WEE = Ntwk.wEE_initial; %rand(Ntwk.Exct.N, Ntwk.Exct.N);
% - spiking events
Espikes = gpuArray.zeros(Ntwk.Exct.N,1);
Ispikes = gpuArray.zeros(Ntwk.Inhbt.N,1);
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
% - intermediate variable for STPD convolution over time
xEpre = gpuArray.zeros(1,Ntwk.Exct.N);
xEpost = gpuArray.zeros(Ntwk.Exct.N,1);
xIpre = gpuArray.zeros(1,Ntwk.Inhbt.N);
xIpost = gpuArray.zeros(Ntwk.Inhbt.N,1);
% Dynamic variables of the example neurons
exampleE = Ntwk.Smpl.E; % E1, E2, EShare1, EShare2
exampleI = Ntwk.Smpl.I; % I1, I2, IShare
Exmpl.ExctV = [ExctV(exampleE,1)'; zeros(timesteps-1, numel(exampleE))];
Exmpl.InhbtV = [InhbtV(exampleI,1)'; zeros(timesteps-1, numel(exampleI))];
Exmpl.ExctgE = [ExctgE(exampleE,1)'; zeros(timesteps-1, numel(exampleE))];
Exmpl.ExctgI = [ExctgI(exampleE,1)'; zeros(timesteps-1, numel(exampleE))];
Exmpl.InhbtgE = [InhbtgE(exampleI,1)'; zeros(timesteps-1, numel(exampleI))];
Exmpl.InhbtgI = [InhbtgI(exampleI,1)'; zeros(timesteps-1, numel(exampleI))];
Exmpl.WEI = gpuArray.zeros(timesteps, numel(exampleI), numel(exampleE));
Exmpl.WEI(1,:,:) = WEI(exampleI,exampleE);
Exmpl.WIE = gpuArray.zeros(timesteps, numel(exampleE), numel(exampleI));
Exmpl.WIE(1,:,:) = WIE(exampleE,exampleI);
Exmpl.WEE = gpuArray.zeros(timesteps, numel(exampleE), numel(exampleE));
Exmpl.WEE(1,:,:) = WEE(exampleE,exampleE);
Exmpl.Espikes = [Espikes(exampleE,1)'; zeros(timesteps-1, numel(exampleE))];
Exmpl.Ispikes = [Ispikes(exampleI,1)'; zeros(timesteps-1, numel(exampleI))];
Exmpl.xEpre = [xEpre(1,exampleE); zeros(timesteps-1, numel(exampleE))];
Exmpl.xEpost = [xEpost(exampleE,1)'; zeros(timesteps-1, numel(exampleE))];
Exmpl.xIpre = [xIpre(1,exampleI); zeros(timesteps-1, numel(exampleI))];
Exmpl.xIpost = [xIpost(exampleI,1)'; zeros(timesteps-1, numel(exampleI))];

% Other parameters
Ib = 150; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
ExctNoise = Ib + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
InhbtNosie = Ib + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
Exmpl.ExctNoise = [ExctNoise(exampleE,1)'; zeros(timesteps-1, numel(exampleE))]; % ExctNoise(exampleE,1)
% simulation start...
filename = 'RealtimeMonitor';
moviename = 'RealtimeMonitor'; % Name of the video file
% Prepare the video file
writerObj = VideoWriter(fullfile(plotdir, moviename));
writerObj.FrameRate = 100; % Adjust frame rate as needed
writerObj.Quality = 100;   % Set quality to maximum for best results (only for MPEG-4)
open(writerObj);
ExctVec = 1:Ntwk.Exct.N;
InhbtVec = [1:Ntwk.Inhbt.N]*4;
h1 = figure;
h2 = figure('Position', [1, 1, 1920, 960]);hold on;
pbaspect([2 1 1]);
h2.Color = 'k'; % Sets the figure background to black
ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
xlabel('x (\mum)', 'Color', 'w');
ylabel('y (\mum)', 'Color', 'w');
xlim([-Ntwk.Scale, Ntwk.Scale]);
ylim([-Ntwk.Scale/2, Ntwk.Scale/2]);
for t = 1:(timesteps-1)
    % input spikes
    InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
    InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(1,t);
    InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(2,t);
    % Synaptic plasticity
    xEpre = xEpre -(xEpre/Ntwk.ExctSTDP.tau_prepost)*dt + Espikes';
    Exmpl.xEpre(t+1,:) = xEpre(exampleE);
    xEpost = xEpost -(xEpost/Ntwk.ExctSTDP.tau_postpre)*dt + Espikes;
    Exmpl.xEpost(t+1,:) = xEpost(exampleE)';
    xIpre = xIpre -(xIpre/Ntwk.InhbtSTDP.tau_prepost)*dt + Ispikes';
    Exmpl.xIpre(t+1,:) = xIpre(exampleI);
    xIpost = xIpost -(xIpost/Ntwk.InhbtSTDP.tau_postpre)*dt + Ispikes;
    Exmpl.xIpost(t+1,:) = xIpost(exampleI)';
    WEI = WEI + Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xIpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Ispikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Ispikes*xEpre + Ntwk.ExctSTDP.intercept_post*Ispikes*ones(size(Espikes'))); % pre -> post
    Exmpl.WEI(t+1,:,:) = WEI(exampleI, exampleE);
    WIE = WIE + Ntwk.InhbtSTDP.eta*(Ntwk.InhbtSTDP.sign_postpre*xEpost*Ispikes' + Ntwk.InhbtSTDP.intercept_pre*ones(size(Espikes))*Ispikes' ... % post -> pre
        + Ntwk.InhbtSTDP.sign_prepost*Espikes*xIpre + Ntwk.InhbtSTDP.intercept_post*Espikes*ones(size(Ispikes'))); % pre -> post
    Exmpl.WIE(t+1,:,:) = WIE(exampleE, exampleI);
    WEE = WEE + Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xEpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Espikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Espikes*xEpre + Ntwk.ExctSTDP.intercept_post*Espikes*ones(size(Espikes'))); % pre -> post
    Exmpl.WEE(t+1,:,:) = WEE(exampleE, exampleE);

    % Synaptic activities
    ExctgE = ExctgE - ExctgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_Input*InputSpikes + Ntwk.Cnnct_EE.*WEE*Espikes);
    Exmpl.ExctgE(t+1,:) = ExctgE(exampleE)';
    ExctgI = ExctgI - ExctgI/Ntwk.Synapse.tauInhbt*dt ... % inhibitory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarI*(Ntwk.Cnnct_IE.*WIE*Ispikes);
    Exmpl.ExctgI(t+1,:) = ExctgI(exampleE)';
    InhbtgE = InhbtgE - InhbtgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Inhbt neurons
        + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_EI.*WEI*Espikes);
    Exmpl.InhbtgE(t+1,:) = InhbtgE(exampleI)';
    InhbtgI = InhbtgI - InhbtgI/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Inhbt neurons
    Exmpl.InhbtgI(t+1,:) = InhbtgI(exampleI)';

    % Updating OU noise
    ExctNoise = ExctNoise + ((Ib - ExctNoise)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm)*dt;
    Exmpl.ExctNoise(t+1,:) = ExctNoise(exampleE)';
    InhbtNosie = InhbtNosie + ((Ib - InhbtNosie)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm)*dt;
    % Membrane potential change for Exct neurons
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + ExctgE.*(Ntwk.VE - ExctV) + ExctgI.*(Ntwk.VI - ExctV) + ExctNoise)/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    ExctRefraction = ExctRefraction - 1;
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    ExctV(Espikes) = Ntwk.Vfire;
    Exmpl.ExctV(t+1,:) = ExctV(exampleE);
    Exmpl.Espikes(t+1,:) = Espikes(exampleE);
    ExctRefraction(Espikes) = refractionPeriod.E;

    % Membrane potential change for Inhbt neurons
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + InhbtgE.*(Ntwk.VE - InhbtV) + InhbtgI.*(Ntwk.VI - InhbtV) + InhbtNosie)/Ntwk.Inhbt.Cm*dt;
    InhbtV = InhbtV + dV;
    InhbtV(InhbtRefraction>0) = Ntwk.Vreset;
    InhbtRefraction = InhbtRefraction - 1;
    % Inhbt neurons fire
    Ispikes = InhbtV > Ntwk.Vth;
    InhbtV(Ispikes) = Ntwk.Vfire-10;
    Exmpl.InhbtV(t+1,:) = InhbtV(exampleI)';
    Exmpl.Ispikes(t+1,:) = Ispikes(exampleI)';
    InhbtRefraction(Ispikes) = refractionPeriod.I;
    % visualization
    figure(h2);
    plot(Ntwk.Exct.Location(Espikes,1), Ntwk.Exct.Location(Espikes,2),'w.', 'MarkerSize', 2);
    plot(Ntwk.Inhbt.Location(Ispikes,1), Ntwk.Inhbt.Location(Ispikes,2),'r.', 'MarkerSize', 2);
    if mod(t*dt, 10) == 1
        drawnow; % Update figure window
        % Capture the plot as an image and write to video
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
        cla;
    end
    if t <= 5000/dt
        figure(h1);
        subplot(1,2,1); hold on;
        plot(time(t)*ones(sum(Espikes),1), ExctVec(Espikes), 'k.', 'MarkerSize', 2);
        plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 2);
    end
    if t >= timesteps-1 - 5000/dt
        figure(h1);
        subplot(1,2,2); hold on;
        plot(time(t)*ones(sum(Espikes),1), ExctVec(Espikes), 'k.', 'MarkerSize', 2);
        plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 2);
    end
    if mod(t*dt, 10) == 0
        fprintf('.');
        if mod(t*dt,100)==0
            fprintf('%1.4f s\n', t*dt/1000);
        end
    end
end
fprintf('Successfully completed\n');
% Close the video file
close(writerObj);
figure(h1); subplot(1,2,1);
axis([leftt rightt 0 Ntwk.Exct.N]);
% axis([43575, 44175, 0 Ntwk.Exct.N]);
xlabel('Time (ms)');
ylabel('Neurons');
title('Raster plot of Exct/Inhbt neurons');
mysavefig(h, filename, plotdir, 12, [8, 6], 1);
subplot(1,2,2);
axis([duration-5000, duration, 0 Ntwk.Exct.N]);
xlabel('Time (ms)');
ylabel('Neurons');
title('Raster plot of Exct/Inhbt neurons');
mysavefig(h, filename, plotdir, 12, [8, 6], 1);
savefig(h, fullfile(plotdir,filename));

%% visualizing example neurons
% load(Rsltfile);
h = figure;
filename = 'Example neurons activity';
subplot(2,1,1); hold on;
lg = [];
lg(1) = plot(time/1000, Exmpl.ExctgE(:,1), 'k-');
lg(2) = plot(time/1000, Exmpl.ExctgI(:,1), 'k--');
lg(3) = plot(time/1000, Exmpl.InhbtgE(:,1), 'r-');
lg(4) = plot(time/1000, Exmpl.InhbtgI(:,1), 'r--');
legend(lg, {'gE on Exct cell','gI on Exct cell', 'gE on Inhbt cell','gI on Inhbt cell'}, 'Location','best');
xlim([0 duration/1000]);
xlabel('Time (s)');
ylabel('Conductance (nS)');
title('Synaptic activity');
mysavefig(h, filename, plotdir, 12, [8,4], 1);
subplot(2,2,3); hold on;
lg = [];
lg(1) = plot(time, Exmpl.ExctgE(:,1), 'k-');
lg(2) = plot(time, Exmpl.ExctgI(:,1), 'k--');
lg(3) = plot(time, Exmpl.InhbtgE(:,1), 'r-');
lg(4) = plot(time, Exmpl.InhbtgI(:,1), 'r--');
xlim([leftt, leftt+700]);
xlabel('Time (ms)');
ylabel('Conductance (nS)');
title('Synaptic activity');
mysavefig(h, filename, plotdir, 12, [8,4], 1);
subplot(2,2,4); hold on;
lg = [];
lg(1) = plot(time, Exmpl.ExctV(:,1), 'k-');
lg(2) = plot(time, Exmpl.InhbtV(:,1), 'r-');
legend(lg, {'Exct', 'Inhbt'}, 'Location','best');
xlim([leftt, leftt+700]);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Activity of a single neuron');
mysavefig(h, filename, plotdir, 12, [8, 4], 1);
%% Synaptic weights on example neurons
h = figure;
filename = 'Example synaptic weights';
subplot(4,2,1); hold on;
lg = [];
lg(1) = plot(time/1000, Exmpl.xEpre(:,1), 'k-');
lg(2) = plot(time/1000, Exmpl.xIpre(:,1), 'r-');
%lg(3) = plot(time/1000, Exmpl.xEpost, 'k--');
%lg(4) = plot(time/1000, Exmpl.xIpost, 'r--');
xlim([0 duration/1000]);
xlabel('Time (s)');
ylabel('Integration');
title('STDP Integration');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,3); hold on;
plot(time/1000, squeeze(Exmpl.WEI(:,1,1)), 'k-');
xlim([0 duration/1000]);
xlabel('Time (s)');
ylabel('wEI');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,5); hold on;
plot(time/1000, squeeze(Exmpl.WIE(:,1,1)), 'r-');
xlim([0 duration/1000]);
xlabel('Time (s)');
ylabel('wIE');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,7); hold on;
plot(time/1000, squeeze(Exmpl.WEE(:,2,1)), 'k--');
xlim([0 duration/1000]);
xlabel('Time (s)');
ylabel('wEE');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,2); hold on;
lg = [];
lg(1) = plot(time, Exmpl.xEpre(:,1), 'k-');
lg(2) = plot(time, Exmpl.xIpre(:,1), 'r-');
% lg(3) = plot(time, Exmpl.xEpost(1,:), 'k--');
% lg(4) = plot(time, Exmpl.xIpost(1,:), 'r--');
legend(lg, {'xE', 'xI'}, 'Location','best');
xlim([leftt leftt+700]);
xlabel('Time (ms)');
ylabel('Integration');
title('STDP Integration');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,4); hold on;
plot(time, squeeze(Exmpl.WEI(:,1,1)), 'k-');
legend('WEI', 'Location','best');
xlim([leftt leftt+700]);
xlabel('Time (ms)');
ylabel('wEI');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,6); hold on;
plot(time, squeeze(Exmpl.WIE(:,1,1)), 'r-');
legend('WIE', 'Location','best');
xlim([leftt leftt+700]);
xlabel('Time (ms)');
ylabel('wIE');
mysavefig(h, filename, plotdir, 12, [8,8], 1);
subplot(4,2,8); hold on;
plot(time, squeeze(Exmpl.WEE(:,2,1)), 'k--');
legend('WEE', 'Location','best');
xlim([leftt leftt+700]);
xlabel('Time (ms)');
ylabel('wEE');
mysavefig(h, filename, plotdir, 12, [8,8], 1);

%% visualize weight change
h = figure;
filename = sprintf('WEI_change_%1.1fs', duration/1000);
imagesc(WEI - Ntwk.wEI_initial)
colormap(bluewhitered);
c = colorbar;
c.Label.String = 'Weight change';
c.Location = 'northoutside';
xlabel("Exct neurons");
ylabel("Inhbt neurons");
mysavefig(h, filename, plotdir, 12, [2.5, 2.81], 1);

h = figure;
filename = sprintf('WIE_change_%1.1fs', duration/1000);
imagesc(WIE - Ntwk.wIE_initial)
colormap(bluewhitered);
c = colorbar;
c.Label.String = 'Weight change';
c.Location = 'northoutside';
xlabel("Inhbt neurons");
ylabel("Exct neurons");
mysavefig(h, filename, plotdir, 12, [2.5, 2.81], 1);

h = figure;
filename = sprintf('WEE_change_%1.1fs', duration/1000);
imagesc(WEE - Ntwk.wEE_initial)
colormap(bluewhitered);
c = colorbar;
c.Label.String = 'Weight change';
c.Location = 'northoutside';
xlabel("Exct neurons");
ylabel("Exct neurons");
mysavefig(h, filename, plotdir, 12, [2.5, 2.81], 1);
%%
EvalTuning(Ntwk,WEE,WEI,WIE,OKeeffe,plotdir);

%% Save results
close all;
clearvars -except 'Ntwk' 'Seq' 'Exmpl' 'WEI' 'WIE' 'WEE' 'Rsltfile';
save(Rsltfile, 'Ntwk', 'Seq', 'Exmpl', 'WEI', 'WIE', 'WEE', '-v7.3');
