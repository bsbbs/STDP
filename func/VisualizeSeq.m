function VisualizeSeq(evs, values, OKeeffe, plotdir)

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
    % Plot the selected input boxcar
    subplot(3,1,2); hold on;
    timevec = 1:(evs(6,1)*1000/dt);
    for ii = Ntwk.Input.Source:-1:1
        plot(timevec*dt/1000, Seq(timevec,ii)*.9+(ii), '-', "Color", OKeeffe(ii,:), 'LineWidth', 1);
    end
    title('Input signal');
    xlabel('Time (s)');
    ylabel('Channel');
    axis([0 evs(6,1), .5 Ntwk.Input.Source*1.45]); % Adjust the axis for better visualization
    % ylim([.5, Ntwk.Input.Source*1.45]);
    yticks([1:Ntwk.Input.Source]);
    mysavefig(h, filename, plotdir, 12, [3,6], 1);
    % Plot the raster plot
    subplot(3,1,3);
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
end