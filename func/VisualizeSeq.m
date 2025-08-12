function VisualizeSeq(Inputstruct, Seq, dt, OKeeffe, plotdir)
%% Visualize the Sgnlsuences of inputs within representative ranges
evs = Seq.evs;
values = Seq.values;
ExmplTrials = min(8, size(evs,1));
time = dt/1000:dt/1000:(max(evs(ExmplTrials,:))+1);
Sgnls = zeros(numel(time),size(evs,2));
for ii = 1:size(evs,2)
    for evi = 1:size(evs,1)
        Sgnls(time >= evs(evi,ii) & time<=evs(evi,ii)+.5,ii) = values(evi,ii);
    end
end
h = figure;
filename = 'InputSequences';
subplot(3,1,1); hold on;
for ii = Inputstruct.Source:-1:1
    plot(time, Sgnls(:,ii)*.9+(ii), '-', "Color", OKeeffe(ii,:), 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (s)');
ylabel('Channel');
% axis([0 duration/1000, .5, Inputstruct.Source*1.45]); % Adjust the axis for better visualization
% ylim([.5, Inputstruct.Source*1.45]);
yticks([1:Inputstruct.Source]);
mysavefig(h, filename, plotdir, 12, [3,6], 1);
% Simulate multiple neurons under driven of the same input signal
% Parameters for Poisson distribution
% Generate Poisson distributed spikes
% The probability of a spike in each time bin
spikeProbability = Inputstruct.spikeRate * dt;
% Poisson generator: generate random numbers and compare to spike probability
leftt = floor(min(find(Sgnls(:,1)>0, 1 ), find(Sgnls(:,2)>0, 1 ))*dt/100)*100;
rightt = ceil((max(find(Sgnls(:,1)>0, 1 ), find(Sgnls(:,2)>0, 1 ))*dt+500)/100)*100;
tmpsteps = (rightt - leftt)/dt;
if gpuDeviceCount > 0
    InputSpikes = gpuArray.rand(tmpsteps, Inputstruct.N);
else
    InputSpikes = rand(tmpsteps, Inputstruct.N);
end
trunck = 50000;
for i = 1:ceil(tmpsteps/trunck)
    timevec = 1+(i-1)*trunck:min(i*trunck, tmpsteps);
    InputSpikes(timevec, Inputstruct.Origins == 1) = InputSpikes(timevec, Inputstruct.Origins == 1) < spikeProbability*Sgnls(timevec + leftt/dt,1);
    InputSpikes(timevec, Inputstruct.Origins == 2) = InputSpikes(timevec, Inputstruct.Origins == 2) < spikeProbability*Sgnls(timevec + leftt/dt,2);
end
% Plot the raster plot
subplot(3,1,2);
hold on;
smpl = min(Inputstruct.N-1,5);
for n = 1:smpl:Inputstruct.N
    % Find indices where spikes occur
    spikeIndices = find(InputSpikes(:, n));
    % Convert indices to time
    spikeTimes = time(spikeIndices+leftt/dt);
    % Plot spikes for this neuron
    if Inputstruct.Origins(n) == 1
        plot(spikeTimes, n * ones(size(spikeTimes)), '.', "Color", OKeeffe(1,:), 'MarkerSize', 2);
    elseif Inputstruct.Origins(n) == 2
        plot(spikeTimes, n * ones(size(spikeTimes)), '.', "Color", OKeeffe(2,:), 'MarkerSize', 2);
    end
end
xlabel('Time (ms)');
ylabel('Input neurons');
title('Spike train');
axis([leftt/1000 rightt/1000 0.5 Inputstruct.N+0.5]); % Adjust the axis for better visualization
mysavefig(h, filename, plotdir, 12, [3,6], 1);
end