function RunnerNMDA1D(Ntwk, Seq, dt, Runningdir)
%% Simulating the neural network
% Initializing the network status at t=0
% - membrane potentials
ExctV = Ntwk.VL*gpuArray.ones(Ntwk.Exct.N,1);
InhbtV = Ntwk.VL*gpuArray.ones(Ntwk.Inhbt.N,1);
ExctRefraction = gpuArray.zeros(Ntwk.Exct.N,1);
InhbtRefraction = gpuArray.zeros(Ntwk.Inhbt.N,1);
% - synaptic conductance (presynaptic-activity dependent)
InputxNMDAi = gpuArray.zeros(Ntwk.Input.N, 1);
InputgNMDAi = gpuArray.zeros(Ntwk.Input.N, 1);
ExctgAMPA = gpuArray.zeros(Ntwk.Exct.N, 1);
ExNMDAi = gpuArray.zeros(Ntwk.Exct.N, 1);
EgNMDAi = gpuArray.zeros(Ntwk.Exct.N, 1);
ExctgNMDA = gpuArray.zeros(Ntwk.Exct.N, 1);
ExctgGABA = gpuArray.zeros(Ntwk.Exct.N, 1);
InhbtgAMPA = gpuArray.zeros(Ntwk.Inhbt.N, 1);
IxNMDAi = gpuArray.zeros(Ntwk.Exct.N, 1);
IgNMDAi = gpuArray.zeros(Ntwk.Exct.N, 1);
InhbtgNMDA = gpuArray.zeros(Ntwk.Inhbt.N, 1);
InhbtgGABA = gpuArray.zeros(Ntwk.Inhbt.N, 1);
% - synaptic weights (plastic according to STDP rules)
WEI = Ntwk.wEI_initial; % rand(Ntwk.Inhbt.N, Ntwk.Exct.N);
WIE = Ntwk.wIE_initial; % rand(Ntwk.Exct.N, Ntwk.Inhbt.N);
WEE = Ntwk.wEE_initial; % rand(Ntwk.Exct.N, Ntwk.Exct.N);
% - spiking events
Espikes = gpuArray.zeros(Ntwk.Exct.N,1);
Ispikes = gpuArray.zeros(Ntwk.Inhbt.N,1);
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
spikeProbability = Ntwk.Input.spikeRate * dt;
% - intermediate variable for STPD convolution over time
xEpre = gpuArray.zeros(1,Ntwk.Exct.N);
xEpost = gpuArray.zeros(Ntwk.Exct.N,1);
xIpre = gpuArray.zeros(1,Ntwk.Inhbt.N);
xIpost = gpuArray.zeros(Ntwk.Inhbt.N,1);
% Other parameters
Ib = 150; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
ExctNoise = Ib + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
InhbtNosie = Ib + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
% % Dynamic variables of the example neurons
% smplintrvl = dt; % ms
smplonsets = round(Seq.evs(round(linspace(1,numel(Seq.evs(:,1)),30)),1)*1000/dt);
smplonsets = unique(smplonsets);
% display(smplonsets*dt/1000);
% smplsteps = (500+200)/smplintrvl;
% timevec = nan(smplsteps-1,1);
% smplE = Ntwk.Smpl.E; % E1, E2, EShare1, EShare2
% smplI = Ntwk.Smpl.I; % I1, I2, IShare
% Smpl.ExctV = [ExctV(smplE,1)'; nan(smplsteps-1, numel(smplE))];
% Smpl.InhbtV = [InhbtV(smplI,1)'; nan(smplsteps-1, numel(smplI))];
% Smpl.ExctgAMPA = [ExctgAMPA(smplE,1)'; nan(smplsteps-1, numel(smplE))];
% Smpl.ExctgNMDA = [ExctgNMDA(smplE,1)'; nan(smplsteps-1, numel(smplE))];
% Smpl.ExctgGABA = [ExctgGABA(smplE,1)'; nan(smplsteps-1, numel(smplE))];
% Smpl.InhbtgAMPA = [InhbtgAMPA(smplI,1)'; nan(smplsteps-1, numel(smplI))];
% Smpl.InhbtgNMDA = [InhbtgNMDA(smplI,1)'; nan(smplsteps-1, numel(smplI))];
% Smpl.InhbtgGABA = [InhbtgGABA(smplI,1)'; nan(smplsteps-1, numel(smplI))];
% Smpl.WEI = gpuArray.zeros(smplsteps, numel(smplI), numel(smplE));
% Smpl.WEI(1,:,:) = WEI(smplI, smplE);
% Smpl.WIE = gpuArray.zeros(smplsteps, numel(smplE), numel(smplI));
% Smpl.WIE(1,:,:) = WIE(smplE, smplI);
% Smpl.WEE = gpuArray.zeros(smplsteps, numel(smplE), numel(smplE));
% Smpl.WEE(1,:,:) = WEE(smplE, smplE);
% Smpl.Espikes = [Espikes(smplE,1)'; zeros(smplsteps-1, numel(smplE))];
% Smpl.Ispikes = [Ispikes(smplI,1)'; zeros(smplsteps-1, numel(smplI))];
% Smpl.xEpre = [xEpre(1,smplE); zeros(smplsteps-1, numel(smplE))];
% Smpl.xEpost = [xEpost(smplE,1)'; zeros(smplsteps-1, numel(smplE))];
% Smpl.xIpre = [xIpre(1,smplI); zeros(smplsteps-1, numel(smplI))];
% Smpl.xIpost = [xIpost(smplI,1)'; zeros(smplsteps-1, numel(smplI))];
% Smpl.ExctNoise = [ExctNoise(smplE,1)'; zeros(smplsteps-1, numel(smplE))]; % ExctNoise(exampleE,1)
% % prepare for visualization
% if visualize == 1
%     h1 = figure; hold on;
%     xlabel('Time (ms)');
%     ylabel('Neurons');
%     ylim([-Ntwk.XScale, Ntwk.XScale]);
%     title('Raster plot of Exct/Inhbt neurons');
%     h2 = figure('Position', [1, 1, 1920, 460]); hold on;
%     pbaspect([2 1 1]);
%     h2.Color = 'k'; % Sets the figure background to black
%     ax = gca;
%     ax.Color = 'k';
%     ax.XColor = 'w';
%     ax.YColor = 'w';
%     xlabel('x (\mum)', 'Color', 'w');
%     ylabel('y (\mum)', 'Color', 'w');
%     xlim([-Ntwk.XScale, Ntwk.XScale]);
%     ylim([-Ntwk.YScale, Ntwk.YScale]);
% end
% Save spike train to incoorperate the delays
bankwidth = round(max([Ntwk.Delay.EE, Ntwk.Delay.EE, Ntwk.Delay.IE])/dt);
SpkTrnE = nan(Ntwk.Exct.N+1,bankwidth); tickE = 0;
SpkTrnI = nan(Ntwk.Inhbt.N+1,bankwidth); tickI = 0;
% simulation start...
% smpl = 0;
duration = max(Seq.evs(:))+1; % secs
timesteps = round(duration*1000/dt);
for t = 1:(timesteps-1)
    % input spikes, Poission process
    Evis = (Seq.evs*1000 - t*dt) >= 0 & (Seq.evs*1000 - t*dt) <= 500;
    InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
    for inputi = 1:Ntwk.Input.Source
        if any(Evis(:,inputi))
            InputSpikes(Ntwk.Input.Origins == inputi) = InputSpikes(Ntwk.Input.Origins == inputi) < spikeProbability*Seq.values(Evis(:,inputi),inputi);
        else
            InputSpikes(Ntwk.Input.Origins == inputi) = 0;
        end
    end
    % Synaptic plasticity
    xEpre = xEpre -(xEpre/Ntwk.ExctSTDP.tau_prepost)*dt + Espikes';
    xEpost = xEpost -(xEpost/Ntwk.ExctSTDP.tau_postpre)*dt + Espikes;
    xIpre = xIpre -(xIpre/Ntwk.InhbtSTDP.tau_prepost)*dt + Ispikes';
    xIpost = xIpost -(xIpost/Ntwk.InhbtSTDP.tau_postpre)*dt + Ispikes;
    WEI = WEI + Ntwk.Cnnct_EI.*(1-WEI).*(Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xIpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Ispikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Ispikes*xEpre + Ntwk.ExctSTDP.intercept_post*Ispikes*ones(size(Espikes')))); % pre -> post
    WEI(WEI<0) = 0;
    WIE = WIE + Ntwk.Cnnct_IE.*(1-WIE).*(Ntwk.InhbtSTDP.eta*(Ntwk.InhbtSTDP.sign_postpre*xEpost*Ispikes' + Ntwk.InhbtSTDP.intercept_pre*ones(size(Espikes))*Ispikes' ... % post -> pre
        + Ntwk.InhbtSTDP.sign_prepost*Espikes*xIpre + Ntwk.InhbtSTDP.intercept_post*Espikes*ones(size(Ispikes')))); % pre -> post
    WIE(WIE<0) = 0;
    WEE = WEE + Ntwk.Cnnct_EE.*(1-WEE).*(Ntwk.ExctSTDP.eta*(Ntwk.ExctSTDP.sign_postpre*xEpost*Espikes' + Ntwk.ExctSTDP.intercept_pre*ones(size(Espikes))*Espikes' ... % post -> pre
        + Ntwk.ExctSTDP.sign_prepost*Espikes*xEpre + Ntwk.ExctSTDP.intercept_post*Espikes*ones(size(Espikes')))); % pre -> post
    WEE(WEE<0) = 0;

    % Synaptic activities
    % AMPA and NMDA on Excitatory neurons
    ExctgAMPA = ExctgAMPA - ExctgAMPA/Ntwk.Synapse.tauExct*dt; % AMPA on Exct neurons
    InputxNMDAi = InputxNMDAi + (-InputxNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Input neurons
    ExNMDAi = ExNMDAi + (-ExNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Exct neurons
    if any(InputSpikes)
        ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*Ntwk.wInput*InputSpikes;
        InputxNMDAi = InputxNMDAi + InputSpikes;
    end
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - round(Ntwk.Delay.EE/dt));
    if ~isempty(DlyEspikes)
        ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*WEE*DlyEspikes;
        ExNMDAi = ExNMDAi + DlyEspikes;
    end
    InputgNMDAi = InputgNMDAi + (-InputgNMDAi/Ntwk.Synapse.NMDA.taudecay + (1-InputgNMDAi).*InputxNMDAi)*dt;
    EgNMDAi = EgNMDAi + (-EgNMDAi/Ntwk.Synapse.NMDA.taudecay + (1-EgNMDAi).*ExNMDAi)*dt;
    ExctgNMDA = Ntwk.Synapse.gbarE*(WEE*EgNMDAi + Ntwk.wInput*InputgNMDAi);

    % GABA on excitatory neurons
    ExctgGABA = ExctgGABA - ExctgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Exct neurons
    DlyIspikes = SpkTrnI(2:end,SpkTrnI(1,:) == t - round(Ntwk.Delay.IE/dt));
    if ~isempty(DlyIspikes)
        ExctgGABA = ExctgGABA + Ntwk.Synapse.gbarI*WIE*DlyIspikes;
    end

    % AMPA and NMDA on inhibitory neurons
    InhbtgAMPA = InhbtgAMPA - InhbtgAMPA/Ntwk.Synapse.tauExct*dt; % AMPA on Inhbt neurons
    IxNMDAi = IxNMDAi + (-IxNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Exct neurons, targetting inhibitory neurons
    DlyEspikes = SpkTrnE(2:end, SpkTrnE(1,:) == t - round(Ntwk.Delay.EI/dt));
    if ~isempty(DlyEspikes)
        InhbtgAMPA = InhbtgAMPA + Ntwk.Synapse.gbarE*(WEI*DlyEspikes);
        IxNMDAi = IxNMDAi + DlyEspikes;
    end
    IgNMDAi = IgNMDAi + (-IgNMDAi/Ntwk.Synapse.NMDA.taudecay + (1-IgNMDAi).*IxNMDAi)*dt;
    InhbtgNMDA = Ntwk.Synapse.gbarE*WEI*IgNMDAi;

    % GABA on inhibitory neurons
    InhbtgGABA = InhbtgGABA - InhbtgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Inhbt neurons

    % Updating OU noise
    ExctNoise = ExctNoise + ((Ib - ExctNoise)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm)*dt;
    InhbtNosie = InhbtNosie + ((Ib - InhbtNosie)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm)*dt;

    % Membrane potential change for Exct neurons
    INMDA = ExctgNMDA./(1+exp(-Ntwk.Synapse.NMDA.a*ExctV/Ntwk.Synapse.NMDA.b)*Ntwk.Synapse.NMDA.Mg2).*(Ntwk.VE - ExctV);
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + ExctgAMPA.*(Ntwk.VE - ExctV) + INMDA + ExctgGABA.*(Ntwk.VI - ExctV) + ExctNoise)/Ntwk.Exct.Cm*dt;
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
    INMDA = InhbtgNMDA./(1+exp(-Ntwk.Synapse.NMDA.a*InhbtV/Ntwk.Synapse.NMDA.b)*Ntwk.Synapse.NMDA.Mg2).*(Ntwk.VE - InhbtV);
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + InhbtgAMPA.*(Ntwk.VE - InhbtV) + INMDA + InhbtgGABA.*(Ntwk.VI - InhbtV) + InhbtNosie)/Ntwk.Inhbt.Cm*dt;
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

%     % Example traces
%     if any(t == smplonsets - 100/dt)
%         Sgi = find(t == smplonsets - 100/dt);
%         filename = sprintf('RealtimeMonitor_Segment%i', Sgi);
%         if visualize
%             % Prepare the video file
%             profiles = VideoWriter.getProfiles();
%             disp(profiles);
%             writerObj1 = VideoWriter(fullfile(subplotdir, filename), 'Motion JPEG AVI');
%             writerObj1.FrameRate = 100; % Adjust frame rate as needed
%             writerObj1.Quality = 95;   % Set quality to maximum for best results (only for MPEG-4)
%             open(writerObj1);
%         end
%         smplt = 0; % reset sampling time label
%         smpl = 1; % turn flag on, start sampling
%     end
%     if smpl && mod(t*dt, smplintrvl) == 0
%         smplt = smplt + 1;
%         Smpl.xEpre(smplt,:) = xEpre(smplE);
%         Smpl.xEpost(smplt,:) = xEpost(smplE)';
%         Smpl.xIpre(smplt,:) = xIpre(smplI);
%         Smpl.xIpost(smplt,:) = xIpost(smplI)';
%         Smpl.WEI(smplt,:,:) = WEI(smplI, smplE);
%         Smpl.WIE(smplt,:,:) = WIE(smplE, smplI);
%         Smpl.WEE(smplt,:,:) = WEE(smplE, smplE);
%         Smpl.ExctgAMPA(smplt,:) = ExctgAMPA(smplE)';
%         Smpl.ExctgNMDA(smplt,:) = ExctgNMDA(smplE)';
%         Smpl.ExctgGABA(smplt,:) = ExctgGABA(smplE)';
%         Smpl.InhbtgAMPA(smplt,:) = InhbtgAMPA(smplI)';
%         Smpl.InhbtgNMDA(smplt,:) = InhbtgNMDA(smplI)';
%         Smpl.InhbtgGABA(smplt,:) = InhbtgGABA(smplI)';
%         Smpl.ExctV(smplt,:) = ExctV(smplE);
%         Smpl.Espikes(smplt,:) = Espikes(smplE);
%         Smpl.InhbtV(smplt,:) = InhbtV(smplI)';
%         Smpl.Ispikes(smplt,:) = Ispikes(smplI)';
%         Smpl.ExctNoise(smplt,:) = ExctNoise(smplE)';
%         timevec(smplt) = t*dt;
%     end
%     % taping movie
%     if visualize && smpl
%         figure(h1);
%         plot(t*dt/1000*ones(sum(Espikes),1), Ntwk.Exct.Location(Espikes,1), 'k.', 'MarkerSize', 2);
%         plot(t*dt/1000*ones(sum(Ispikes),1), Ntwk.Inhbt.Location(Ispikes,1), 'r.', 'MarkerSize', 3);
%         figure(h2);
%         plot(Ntwk.Exct.Location(Espikes,1), Ntwk.Exct.Location(Espikes,2),'w.', 'MarkerSize', 8);
%         plot(Ntwk.Inhbt.Location(Ispikes,1), Ntwk.Inhbt.Location(Ispikes,2),'r.', 'MarkerSize', 10);
%         if mod(t*dt, 1) == 0
%             drawnow; % Update figure window
%             frame = getframe(gcf);
%             writeVideo(writerObj1, frame);
%             cla;
%         end
%     end
    if any(t == smplonsets + 600/dt)
        s = struct("WEE", WEE, "WEI", WEI, "WIE", WIE);
        Sgi = find(t == smplonsets + 600/dt);
        filename = sprintf('Weights_Seg%i', Sgi);
        % s = struct("Smpl", Smpl, "timevec", timevec, "WEE", WEE, "WEI", WEI, "WIE", WIE);
        save(fullfile(Runningdir, [filename, '.mat']), '-fromstruct', s);
%         if visualize
%             figure(h1);
%             axis([(smplonsets(Sgi)*dt+[-100, 600])/1000 0 Ntwk.Exct.N]);
%             mysavefig(h1, filename, Runningdir, 12, [4, 4], 1);
%             savefig(h1, fullfile(Runningdir,filename));
%             close(h1);
%             cla;
%             figure(h2);
%             cla;
%             close(writerObj1);
%             close(h2);
%         end
%         smpl = 0; % sampling stopped
    end
    if mod(t*dt, 10) == 0
        fprintf('.');
        if mod(t*dt,100)==0
            fprintf('%1.4f s\n', t*dt/1000);
        end
    end
end

% %% visualize weight change
% h = figure;
% filename = sprintf('WEI_change_%1.1fs', duration/1000);
% imagesc(WEI - Ntwk.wEI_initial)
% colormap(bluewhitered);
% c = colorbar;
% c.Label.String = 'Weight change';
% c.Location = 'northoutside';
% xlabel("Exct neurons");
% ylabel("Inhbt neurons");
% mysavefig(h, filename, Runningdir, 12, [2.5, 2.81], 1);
% 
% h = figure;
% filename = sprintf('WIE_change_%1.1fs', duration/1000);
% imagesc(WIE - Ntwk.wIE_initial)
% colormap(bluewhitered);
% c = colorbar;
% c.Label.String = 'Weight change';
% c.Location = 'northoutside';
% xlabel("Inhbt neurons");
% ylabel("Exct neurons");
% mysavefig(h, filename, Runningdir, 12, [2.5, 2.81], 1);
% 
% h = figure;
% filename = sprintf('WEE_change_%1.1fs', duration/1000);
% imagesc(WEE - Ntwk.wEE_initial)
% colormap(bluewhitered);
% c = colorbar;
% c.Label.String = 'Weight change';
% c.Location = 'northoutside';
% xlabel("Exct neurons");
% ylabel("Exct neurons");
% mysavefig(h, filename, Runningdir, 12, [2.5, 2.81], 1);
% 
% %% Overall tuning
% EvalTuning1D(Ntwk,WEE,WEI,WIE,OKeeffe,Runningdir);
% end
