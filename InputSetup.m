%% Inputs of single source, assume this long-range projection only intervenes Exct neurons
%% Define input connection matrix
Ntwk.Input.Source = 2; % Number of input source(s)
Ntwk.Input.N = 1000*Ntwk.Input.Source; % number of the input (excitarory) neurons
Ntwk.Input.Tube = 50; % um, the Gaussuan standard deviation of long-range projection area. All input neurons locate in the circle area
Ntwk.AxonRange.Input = 70; % 130; % um, for each input neuron, the axon connection to the neighbouring E neuron is assumed the same as Gaussian distance delay as E to E
Ntwk.CnnctProb.Input = 0.1; % The connection probability from input neurons are set as the same as E to E
gpurng(2024);
r = gpuArray.randn(Ntwk.Input.N,1)*Ntwk.Input.Tube; % radius of the input fiber distribution
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
p_Input = Ntwk.CnnctProb.Input*exp(-.5*(DstcInput/Ntwk.AxonRange.Input).^2); % probability of physical connection based on distance
gpurng(2026);
Ntwk.Cnnct_Input = p_Input >= gpuArray.rand(size(p_Input)); % projections from input, 0 or 1
Ntwk.wInput = gpuArray.rand(size(Ntwk.Cnnct_Input)).*Ntwk.Cnnct_Input;
clear XInput YInput XE YE DstcInput p_Input;
if show
    Ntwk = InputTuning(Ntwk, OKeeffe, plotdir);
end
%% Visualize input sequences
h = figure;
filename = 'InputSequences';
subplot(3,1,1); hold on;
for ii = Ntwk.Input.Source:-1:1
    plot(time/1000, Seq(:,ii)*.9+(ii), '-', "Color", OKeeffe(ii,:), 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
axis([0 duration/1000 .5, Ntwk.Input.Source*1.45]); % Adjust the axis for better visualization
% ylim([.5, Ntwk.Input.Source*1.45]);
yticks([1:Ntwk.Input.Source]);
mysavefig(h, filename, plotdir, 12, [3,6], 1);
% Simulate multiple neurons under driven of the same input signal
% Parameters for Poisson distribution
spikeRate = 100/1000; % spikes per milisecond
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = spikeRate * dt;
Ntwk.Input.spikeProbability = spikeProbability;
% Poisson generator: generate random numbers and compare to spike probability
leftt = floor(min(find(Seq(:,1)>0, 1 ), find(Seq(:,2)>0, 1 )*dt)/100)*100;
rightt = ceil((max(find(Seq(:,1)>0, 1 ), find(Seq(:,2)>0, 1 ))*dt+500)/100)*100;
tmpsteps = (rightt - leftt)/dt;
InputSpikes = gpuArray.rand(tmpsteps, Ntwk.Input.N);
trunck = 50000;
for i = 1:ceil(tmpsteps/trunck)
    timevec = 1+(i-1)*trunck:min(i*trunck, tmpsteps);
    InputSpikes(timevec, Ntwk.Input.Origins == 1) = InputSpikes(timevec, Ntwk.Input.Origins == 1) < spikeProbability*Seq(timevec + leftt/dt,1);
    InputSpikes(timevec, Ntwk.Input.Origins == 2) = InputSpikes(timevec, Ntwk.Input.Origins == 2) < spikeProbability*Seq(timevec + leftt/dt,2);
end

% Plot the raster plot
subplot(3,1,2);
hold on;
for n = 1:5:Ntwk.Input.N
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(:, n));
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
mysavefig(h, filename, plotdir, 12, [3,6], 1);
clear InputSpikes spikeTimes spikeIndices SumStrength;
save(fullfile(plotdir, 'Ntwk.mat'), "Ntwk");