% Runner

%% Setup for visualization etc
Setup;

%% Build the neural network
show = 1;
NetworkgeneratorPeriodicGPU;

%% Specify project name and output

Rsltfile = fullfile(plotdir,'Rslts.mat');

%% Build the Input structure
InputSetup;

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
WEI = Ntwk.wEI_initial; % rand(Ntwk.Inhbt.N, Ntwk.Exct.N);
WIE = Ntwk.wIE_initial; % rand(Ntwk.Exct.N, Ntwk.Inhbt.N);
WEE = Ntwk.wEE_initial; % rand(Ntwk.Exct.N, Ntwk.Exct.N);
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
% Other parameters
Ib = 150; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
ExctNoise = Ib + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
InhbtNosie = Ib + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
% Dynamic variables of the example neurons
smplintrvl = .1; % ms
smplonsets = round(evs(round(linspace(1,numel(evs(:,1)),9)),1)*1000/dt);
smplsteps = (500+200)/smplintrvl;
smplE = Ntwk.Smpl.E; % E1, E2, EShare1, EShare2
smplI = Ntwk.Smpl.I; % I1, I2, IShare
Smpl.ExctV = [ExctV(smplE,1)'; nan(smplsteps-1, numel(smplE))];
Smpl.InhbtV = [InhbtV(smplI,1)'; nan(smplsteps-1, numel(smplI))];
Smpl.ExctgE = [ExctgE(smplE,1)'; nan(smplsteps-1, numel(smplE))];
Smpl.ExctgI = [ExctgI(smplE,1)'; nan(smplsteps-1, numel(smplE))];
Smpl.InhbtgE = [InhbtgE(smplI,1)'; nan(smplsteps-1, numel(smplI))];
Smpl.InhbtgI = [InhbtgI(smplI,1)'; nan(smplsteps-1, numel(smplI))];
Smpl.WEI = gpuArray.zeros(smplsteps, numel(smplI), numel(smplE));
Smpl.WEI(1,:,:) = WEI(smplI, smplE);
Smpl.WIE = gpuArray.zeros(smplsteps, numel(smplE), numel(smplI));
Smpl.WIE(1,:,:) = WIE(smplE, smplI);
Smpl.WEE = gpuArray.zeros(smplsteps, numel(smplE), numel(smplE));
Smpl.WEE(1,:,:) = WEE(smplE, smplE);
Smpl.Espikes = [Espikes(smplE,1)'; zeros(smplsteps-1, numel(smplE))];
Smpl.Ispikes = [Ispikes(smplI,1)'; zeros(smplsteps-1, numel(smplI))];
Smpl.xEpre = [xEpre(1,smplE); zeros(smplsteps-1, numel(smplE))];
Smpl.xEpost = [xEpost(smplE,1)'; zeros(smplsteps-1, numel(smplE))];
Smpl.xIpre = [xIpre(1,smplI); zeros(smplsteps-1, numel(smplI))];
Smpl.xIpost = [xIpost(smplI,1)'; zeros(smplsteps-1, numel(smplI))];
Smpl.ExctNoise = [ExctNoise(smplE,1)'; zeros(smplsteps-1, numel(smplE))]; % ExctNoise(exampleE,1)
% simulation start...
% prepare the visualization
ExctVec = 1:Ntwk.Exct.N;
InhbtVec = [1:Ntwk.Inhbt.N]*4;
h1 = figure; hold on;
xlabel('Time (ms)');
ylabel('Neurons');
ylim([1, Ntwk.Exct.N]);
title('Raster plot of Exct/Inhbt neurons');
h2 = figure('Position', [1, 1, 1920, 960]); hold on;
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
% Save spike train to incoorperate the delays
bankwidth = max([Ntwk.Delay.EE, Ntwk.Delay.EE, Ntwk.Delay.IE])/dt;
SpkTrnE = nan(Ntwk.Exct.N+1,bankwidth); tickE = 0;
SpkTrnI = nan(Ntwk.Inhbt.N+1,bankwidth); tickI = 0;
smpl = 0;
for t = 1:(timesteps-1)
    % input spikes
    InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
    InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(t,1);
    InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(t,2);
    
    % Synaptic plasticity
    xEpre = xEpre -(xEpre/Ntwk.ExctSTDP.tau_prepost)*dt + Espikes';
    xEpost = xEpost -(xEpost/Ntwk.ExctSTDP.tau_postpre)*dt + Espikes;
    xIpre = xIpre -(xIpre/Ntwk.InhbtSTDP.tau_prepost)*dt + Ispikes';
    xIpost = xIpost -(xIpost/Ntwk.InhbtSTDP.tau_postpre)*dt + Ispikes;
    WEI = WEI + (1-WEI).*(Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xIpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Ispikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Ispikes*xEpre + Ntwk.ExctSTDP.intercept_post*Ispikes*ones(size(Espikes')))); % pre -> post
    WEI(WEI<0) = 0;
    WIE = WIE + (1-WIE).*(Ntwk.InhbtSTDP.eta*(Ntwk.InhbtSTDP.sign_postpre*xEpost*Ispikes' + Ntwk.InhbtSTDP.intercept_pre*ones(size(Espikes))*Ispikes' ... % post -> pre
        + Ntwk.InhbtSTDP.sign_prepost*Espikes*xIpre + Ntwk.InhbtSTDP.intercept_post*Espikes*ones(size(Ispikes')))); % pre -> post
    WIE(WIE<0) = 0;
    WEE = WEE + (1-WEE).*(Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xEpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Espikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Espikes*xEpre + Ntwk.ExctSTDP.intercept_post*Espikes*ones(size(Espikes')))); % pre -> post
    WEE(WEE<0) = 0;
    % Synaptic activities
    ExctgE = ExctgE - ExctgE/Ntwk.Synapse.tauExct*dt; % excitatory synaptic conductance on Exct neurons
    if any(InputSpikes)
        ExctgE = ExctgE + Ntwk.Synapse.gbarE*Ntwk.Cnnct_Input*InputSpikes;
    end
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EE/dt);
    if ~isempty(DlyEspikes)
        ExctgE = ExctgE + Ntwk.Synapse.gbarE*Ntwk.Cnnct_EE.*WEE*DlyEspikes;
    end
    ExctgI = ExctgI - ExctgI/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Exct neurons
    DlyIspikes = SpkTrnI(2:end,SpkTrnI(1,:) == t - Ntwk.Delay.IE/dt);
    if ~isempty(DlyIspikes)
        ExctgI = ExctgI + Ntwk.Synapse.gbarI*Ntwk.Cnnct_IE.*WIE*DlyIspikes;
    end
    InhbtgE = InhbtgE - InhbtgE/Ntwk.Synapse.tauExct*dt; % excitatory synaptic conductance on Inhbt neurons
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EE/dt);
    if ~isempty(DlyEspikes)
        InhbtgE = InhbtgE + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_EI.*WEI*DlyEspikes);
    end
    InhbtgI = InhbtgI - InhbtgI/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Inhbt neurons

    % Updating OU noise
    ExctNoise = ExctNoise + ((Ib - ExctNoise)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm)*dt;
    InhbtNosie = InhbtNosie + ((Ib - InhbtNosie)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm)*dt;
    
    % Membrane potential change for Exct neurons
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + ExctgE.*(Ntwk.VE - ExctV) + ExctgI.*(Ntwk.VI - ExctV) + ExctNoise)/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    ExctRefraction = ExctRefraction - 1;
    
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    if any(Espikes)
        tickE = tickE + 1;
        if tickE > bankwidth
            tickE = 1;
        end
        SpkTrnE(:, tickE) = [t; Espikes];
    end
    ExctV(Espikes) = Ntwk.Vfire;
    ExctRefraction(Espikes) = refractionPeriod.E;

    % Membrane potential change for Inhbt neurons
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + InhbtgE.*(Ntwk.VE - InhbtV) + InhbtgI.*(Ntwk.VI - InhbtV) + InhbtNosie)/Ntwk.Inhbt.Cm*dt;
    InhbtV = InhbtV + dV;
    InhbtV(InhbtRefraction>0) = Ntwk.Vreset;
    InhbtRefraction = InhbtRefraction - 1;
    % Inhbt neurons fire
    Ispikes = InhbtV > Ntwk.Vth;
    if any(Ispikes)
        tickI = tickI + 1;
        if tickI > bankwidth
            tickI = 1;
        end
        SpkTrnI(:, tickI) = [t; Ispikes];
    end
    InhbtV(Ispikes) = Ntwk.Vfire;
    InhbtRefraction(Ispikes) = refractionPeriod.I;

    % Example traces
    if any(t == smplonsets - 100/dt)
        evi = find(t == smplonsets - 100/dt);
        filename = sprintf('RealtimeMonitor_Event%i', evi);
        % Prepare the video file
        profiles = VideoWriter.getProfiles();
        disp(profiles);
        writerObj1 = VideoWriter(fullfile(plotdir, filename), 'Motion JPEG AVI');
        writerObj1.FrameRate = 100; % Adjust frame rate as needed
        writerObj1.Quality = 95;   % Set quality to maximum for best results (only for MPEG-4)
        open(writerObj1);
        timevec = [];
        smplt = 0; % reset sampling time label
        smpl = 1; % turn flag on, start sampling
    end
    if smpl && mod(t*dt, smplintrvl) == 0
        smplt = smplt + 1;
        Smpl.xEpre(smplt,:) = xEpre(smplE);
        Smpl.xEpost(smplt,:) = xEpost(smplE)';
        Smpl.xIpre(smplt,:) = xIpre(smplI);
        Smpl.xIpost(smplt,:) = xIpost(smplI)';
        Smpl.WEI(smplt,:,:) = WEI(smplI, smplE);
        Smpl.WIE(smplt,:,:) = WIE(smplE, smplI);
        Smpl.WEE(smplt,:,:) = WEE(smplE, smplE);
        Smpl.ExctgE(smplt,:) = ExctgE(smplE)';
        Smpl.ExctgI(smplt,:) = ExctgI(smplE)';
        Smpl.InhbtgE(smplt,:) = InhbtgE(smplI)';
        Smpl.InhbtgI(smplt,:) = InhbtgI(smplI)';
        Smpl.ExctV(smplt,:) = ExctV(smplE);
        Smpl.Espikes(smplt,:) = Espikes(smplE);
        Smpl.InhbtV(smplt,:) = InhbtV(smplI)';
        Smpl.Ispikes(smplt,:) = Ispikes(smplI)';
        Smpl.ExctNoise(smplt,:) = ExctNoise(smplE)';
        timevec = [timevec; t*dt];
    end
    % taping movie
    if smpl
        figure(h1);
        plot(time(t)*ones(sum(Espikes),1), ExctVec(Espikes), 'k.', 'MarkerSize', 2);
        plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 3);
        figure(h2);
        plot(Ntwk.Exct.Location(Espikes,1), Ntwk.Exct.Location(Espikes,2),'w.', 'MarkerSize', 2);
        plot(Ntwk.Inhbt.Location(Ispikes,1), Ntwk.Inhbt.Location(Ispikes,2),'r.', 'MarkerSize', 3);
        if mod(t*dt, 1) == 0
            drawnow; % Update figure window
            frame = getframe(gcf);
            writeVideo(writerObj1, frame);
            cla;
        end
    end
    if any(t == smplonsets + 600/dt)
        save(fullfile(plotdir, [filename, '.mat']), 'Smpl','timevec','WEE','WEI','WIE', '-v7.3');
        figure(h1);
        axis([time(smplonsets(evi))+[-100, 600] 0 Ntwk.Exct.N]);
        mysavefig(h1, filename, plotdir, 12, [4, 4], 1);
        savefig(h1, fullfile(plotdir,filename));
        cla;
        figure(h2);
        cla;
        close(writerObj1);
        smpl = 0; % sampling stopped
    end
    if mod(t*dt, 10) == 0
        fprintf('.');
        if mod(t*dt,100)==0
            fprintf('%1.4f s\n', t*dt/1000);
        end
    end
end
fprintf('Successfully completed\n');

%% visualizing example neurons
% % load(Rsltfile);
% h = figure;
% filename = 'Example neurons activity';
% subplot(2,1,1); hold on;
% lg = [];
% lg(1) = plot(time/1000, Smpl.ExctgE(:,1), 'k-');
% lg(2) = plot(time/1000, Smpl.ExctgI(:,1), 'k--');
% lg(3) = plot(time/1000, Smpl.InhbtgE(:,1), 'r-');
% lg(4) = plot(time/1000, Smpl.InhbtgI(:,1), 'r--');
% legend(lg, {'gE on Exct cell','gI on Exct cell', 'gE on Inhbt cell','gI on Inhbt cell'}, 'Location','best');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('Conductance (nS)');
% title('Synaptic activity');
% mysavefig(h, filename, plotdir, 12, [8,4], 1);
% subplot(2,2,3); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.ExctgE(:,1), 'k-');
% lg(2) = plot(time, Smpl.ExctgI(:,1), 'k--');
% lg(3) = plot(time, Smpl.InhbtgE(:,1), 'r-');
% lg(4) = plot(time, Smpl.InhbtgI(:,1), 'r--');
% xlim([leftt, leftt+700]);
% xlabel('Time (ms)');
% ylabel('Conductance (nS)');
% title('Synaptic activity');
% mysavefig(h, filename, plotdir, 12, [8,4], 1);
% subplot(2,2,4); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.ExctV(:,1), 'k-');
% lg(2) = plot(time, Smpl.InhbtV(:,1), 'r-');
% legend(lg, {'Exct', 'Inhbt'}, 'Location','best');
% xlim([leftt, leftt+700]);
% xlabel('Time (ms)');
% ylabel('Membrane potential (mV)');
% title('Activity of a single neuron');
% mysavefig(h, filename, plotdir, 12, [8, 4], 1);
% %% Synaptic weights on example neurons
% h = figure;
% filename = 'Example synaptic weights';
% subplot(4,2,1); hold on;
% lg = [];
% lg(1) = plot(time/1000, Smpl.xEpre(:,1), 'k-');
% lg(2) = plot(time/1000, Smpl.xIpre(:,1), 'r-');
% %lg(3) = plot(time/1000, Exmpl.xEpost, 'k--');
% %lg(4) = plot(time/1000, Exmpl.xIpost, 'r--');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('Integration');
% title('STDP Integration');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,3); hold on;
% plot(time/1000, squeeze(Smpl.WEI(:,1,1)), 'k-');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wEI');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,5); hold on;
% plot(time/1000, squeeze(Smpl.WIE(:,1,1)), 'r-');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wIE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,7); hold on;
% plot(time/1000, squeeze(Smpl.WEE(:,2,1)), 'k--');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wEE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,2); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.xEpre(:,1), 'k-');
% lg(2) = plot(time, Smpl.xIpre(:,1), 'r-');
% % lg(3) = plot(time, Exmpl.xEpost(1,:), 'k--');
% % lg(4) = plot(time, Exmpl.xIpost(1,:), 'r--');
% legend(lg, {'xE', 'xI'}, 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('Integration');
% title('STDP Integration');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,4); hold on;
% plot(time, squeeze(Smpl.WEI(:,1,1)), 'k-');
% legend('WEI', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wEI');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,6); hold on;
% plot(time, squeeze(Smpl.WIE(:,1,1)), 'r-');
% legend('WIE', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wIE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,8); hold on;
% plot(time, squeeze(Smpl.WEE(:,2,1)), 'k--');
% legend('WEE', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wEE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);

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
% SynpseDynamics;
EvalTuning(Ntwk,WEE,WEI,WIE,OKeeffe,plotdir);
%% Save results
close all;
clearvars -except 'Ntwk' 'Seq' 'Exmpl' 'WEI' 'WIE' 'WEE' 'Rsltfile';
save(Rsltfile, 'Ntwk', 'Seq', 'Smpl', 'WEI', 'WIE', 'WEE', '-v7.3');

