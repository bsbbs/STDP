%% Evaluation of the computation of the circuit
%% Define the input and output directories
% DefineIO1D;
% dt = .1;
% Ntrial = 1600;
% ProjectName = sprintf('SyncNMDA1D100_%i', Ntrial);
% plotdir = fullfile(Projdir, ProjectName);
%% 
evi = 18; % intersectional check based on training duration from 0 to 18 period.
plotdirintersect = fullfile(plotdir, sprintf('Intersect%i', evi));
if ~exist(plotdirintersect,'dir')
    mkdir(plotdirintersect);
end
%% load Ntwk file and trained weights
Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
load(Ntwkfile);
if evi >=1
    filename = sprintf('RealtimeMonitor_Event%i', evi);
    load(fullfile(plotdir, [filename, '.mat']));
else
    WEI = Ntwk.wEI_initial;
    WIE = Ntwk.wIE_initial;
    WEE = Ntwk.wEE_initial;
end
%% Visualize connection and pools of samples
Sum1 = sum(Ntwk.wInput(:, Ntwk.Input.Origins == 1), 2);
Sum2 = sum(Ntwk.wInput(:, Ntwk.Input.Origins == 2), 2);
Cnnct1 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 1), 2);
Cnnct2 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 2), 2);
h = figure;
filename = 'Inputtuning_individualE';
plot(Sum1, Sum2, '.', 'MarkerSize', 2);
xlabel('Tuning weight to Input 1');
ylabel('Tuning weight to Input 2');
mysavefig(h, filename, plotdirintersect, 12, [2.5, 2.2161], 1);

Mat = [sum(~Cnnct1 & ~Cnnct2), sum(Cnnct1 & ~Cnnct2);
    sum(Cnnct2 & ~Cnnct1), sum(Cnnct1 & Cnnct2)];
h = figure;
filename = 'InputTuningNneurons';
imagesc(Mat);
set(gca, 'YDir', 'normal');
colormap('gray');
for i = 1:2
    for j = 1:2
    text(i,j, sprintf('%i',Mat(j,i)),'Color','r');
    end
end
ylabel('Tuning to Input 2');
xlabel('Tuning to Input 1');
mysavefig(h, filename, plotdirintersect, 12, [2.5, 2.2161]); 

h = figure; hold on;
filename = 'InputTuning2D';
scatter(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 3, Sum1, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 3, -Sum2, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Tuning weights';
c.Location = 'northoutside';
xlabel('x (\mum)');
ylabel('y (\mum)');
mysavefig(h, filename, plotdirintersect, 12, [5, 3], 0);

h = figure; hold on;
filename = 'InputTuning2DSlcted';
TuneMask1 = Cnnct1 & ~Cnnct2;
TuneMask2 = ~Cnnct1 & Cnnct2;
scatter(Ntwk.Exct.Location(TuneMask1,1), Ntwk.Exct.Location(TuneMask1,2), 3, Sum1(TuneMask1), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(Ntwk.Exct.Location(TuneMask2,1), Ntwk.Exct.Location(TuneMask2,2), 3, -Sum2(TuneMask2), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Tuning weights';
c.Location = 'northoutside';
xlabel('x (\mum)');
ylabel('y (\mum)');
mysavefig(h, filename, plotdirintersect, 12, [5, 3], 0);
%% define the input values and sequence for testing purpose only
% Seqfile = fullfile(plotdir, 'SeqTest.mat');
% TestN = 16;
% if ~exist(Seqfile, 'file')
%     sgm = 0;
%     values = ParetoSequence(TestN, sgm);
%     [Seq, evs] = Generator(TestN, values, dt);
%     Seq = Seq(:, 1:2);
%     % Time vector
%     duration = length(Seq)*dt; % ms
%     time = [dt:dt:duration]';
%     timesteps = numel(time);
%     save(Seqfile, 'Seq','evs','values','time','duration','timesteps');
% else
%     load(Seqfile);
% end
% h = figure;
% filename = 'TestInputDynamic';
% hold on;
% for ii = 2:-1:1
%     plot(time/1000, Seq(:,ii)*.9+(ii), '-', 'LineWidth', 1);
% end
% title('Input signal');
% xlabel('Time (s)');
% ylabel('Channel');
% axis([0, time(end)/1000, .5, 2*1.45]); % Adjust the axis for better visualization
% yticks([1:2]);
% mysavefig(h, filename, plotdir, 12, [4,2], .1);

duration = 1000; % ms
time = [dt:dt:duration]';
timesteps = numel(time);
% test the network
V1vec = .5:.5:2;
V2vec = flip(1./linspace(1/3,1,11)-1);
Modulation = nan(numel(V1vec),numel(V2vec), 2);
for v1i = 1:numel(V1vec)
    for v2i = 1:numel(V2vec)
        V1 = V1vec(v1i);
        V2 = V2vec(v2i);
        filename = sprintf('RealtimeMonitor_Event%i_testV1%1.1f_V2%1.1f', evi, V1, V2);
        if ~exist(fullfile(plotdirintersect, [filename, '.mat']), 'file')
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
            % - spiking events
            Espikes = gpuArray.zeros(Ntwk.Exct.N,1);
            Ispikes = gpuArray.zeros(Ntwk.Inhbt.N,1);
            refractionPeriod.E = Ntwk.Exct.tauREF/dt;
            refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
            % Other parameters
            Ib = 150; % Vogels et al., 2011; 129; % pA, baseline input to every excitatory neurons
            ExctNoise = Ib + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm; % OU noise on excitatory neurons
            InhbtNosie = Ib + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm; % OU noise on inhibitory neurons
            % Save spike train to incoorperate the delays
            % bankwidth = max([Ntwk.Delay.EE, Ntwk.Delay.EE, Ntwk.Delay.IE])/dt;
            bankwidth = round(max([Ntwk.Delay.EE, Ntwk.Delay.EE, Ntwk.Delay.IE])/dt);
            SpkTrnE = nan(Ntwk.Exct.N+1,bankwidth); tickE = 0;
            SpkTrnI = nan(Ntwk.Inhbt.N+1,bankwidth); tickI = 0;
            h1 = figure; hold on;
            xlabel('Time (ms)');
            ylabel('Neurons');
            ylim([1, Ntwk.Exct.N]);
            title('Raster plot of Exct/Inhbt neurons');
            ExctVec = 1:Ntwk.Exct.N;
            InhbtVec = [1:Ntwk.Inhbt.N]*4;
            % Output variables
            MeanFR = [];
            timevec = [];
            for t = 1:(timesteps-1)
                % input spikes
                InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
                % InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(t,1);
                % InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(t,2);
                if t*dt <= 1000
                    InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < Ntwk.Input.spikeProbability*V1;
                    InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < Ntwk.Input.spikeProbability*V2;
                else
                    InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < Ntwk.Input.spikeProbability*0;
                    InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < Ntwk.Input.spikeProbability*0;
                end
                % Synaptic activities
                % AMPA and NMDA on Excitatory neurons
                ExctgAMPA = ExctgAMPA - ExctgAMPA/Ntwk.Synapse.tauExct*dt; % AMPA on Exct neurons
                InputxNMDAi = InputxNMDAi + (-InputxNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Input neurons
                ExNMDAi = ExNMDAi + (-ExNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Exct neurons
                if any(InputSpikes)
                    ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*Ntwk.wInput*InputSpikes;
                    InputxNMDAi = InputxNMDAi + InputSpikes;
                end
                DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EE/dt);
                if ~isempty(DlyEspikes)
                    ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*WEE*DlyEspikes;
                    ExNMDAi = ExNMDAi + DlyEspikes;
                end
                InputgNMDAi = InputgNMDAi + (-InputgNMDAi/Ntwk.Synapse.NMDA.taudecay + (1-InputgNMDAi).*InputxNMDAi)*dt;
                EgNMDAi = EgNMDAi + (-EgNMDAi/Ntwk.Synapse.NMDA.taudecay + (1-EgNMDAi).*ExNMDAi)*dt;
                ExctgNMDA = Ntwk.Synapse.gbarE*(WEE*EgNMDAi + Ntwk.wInput*InputgNMDAi);

                % GABA on excitatory neurons
                ExctgGABA = ExctgGABA - ExctgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Exct neurons
                DlyIspikes = SpkTrnI(2:end,SpkTrnI(1,:) == t - Ntwk.Delay.IE/dt);
                if ~isempty(DlyIspikes)
                    ExctgGABA = ExctgGABA + Ntwk.Synapse.gbarI*WIE*DlyIspikes;
                end

                % AMPA and NMDA on inhibitory neurons
                InhbtgAMPA = InhbtgAMPA - InhbtgAMPA/Ntwk.Synapse.tauExct*dt; % AMPA on Inhbt neurons
                IxNMDAi = IxNMDAi + (-IxNMDAi/Ntwk.Synapse.NMDA.taurise)*dt; % NMDA from Exct neurons, targetting inhibitory neurons
                DlyEspikes = SpkTrnE(2:end, SpkTrnE(1,:) == t - Ntwk.Delay.EI/dt);
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

                plot(time(t)*ones(sum(Espikes),1), ExctVec(Espikes), 'k.', 'MarkerSize', 2);
                plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 3);
                % Example traces
                MeanFR(t, 1) = sum(Espikes(TuneMask1))/sum(TuneMask1);
                MeanFR(t, 2) = sum(Espikes(TuneMask2))/sum(TuneMask2);
                timevec = [timevec; t*dt];
                if mod(t*dt, 10) == 0
                    fprintf('.');
                    if mod(t*dt,100)==0
                        fprintf('%1.4f s\n', t*dt/1000);
                    end
                end
            end
            fprintf('Successfully completed\n');
            mysavefig(h1, filename, plotdirintersect, 12, [4, 4], 1);
            savefig(h1, fullfile(plotdirintersect,[filename, '.fig']));
            save(fullfile(plotdirintersect, [filename, '.mat']), 'timevec','MeanFR');
        else
            load(fullfile(plotdirintersect, [filename, '.mat']));
        end
        %
        h = figure; hold on;
        filename = sprintf('MeanFRDynmc_V1%1.0fV2%.1f', V1, V2);
        [FR, timep] = PSTH(MeanFR, dt);
        plot(timep, FR(:,1),'-', 'Color', OKeeffe(1,:), 'LineWidth', 2);
        plot(timep, FR(:,2),'-', 'Color', OKeeffe(2,:), 'LineWidth', 2);
        xlabel('Time (ms)');
        ylabel('Firing rate (Hz)');
        mysavefig(h, filename, plotdirintersect, 12, [4, 2], 1);
        Modulation(v1i,v2i, :) = mean(FR(timep>200 & timep<1000,:),1);
    end
end
h = figure; 
filename = 'ModulationEffect';
mygray = gray(numel(V1vec)+1);
subplot(1,2,1);hold on; % direct effect
for v1i = 1:numel(V1vec)
    plot(V2vec,Modulation(v1i,:,2), 'k.', 'Color', mygray(v1i,:), 'MarkerSize', 8);
end
legend({'V_1 = 0.5','V_1 = 1.0','V_1 = 1.5','V_1 = 2.0'},'Location','best');
ylabel('Mean firing rate (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdirintersect, 12, [6, 2.5], 1);
subplot(1,2,2);hold on; % context effect
for v1i = 1:numel(V1vec)
    plot(V2vec,Modulation(v1i,:,1), 'k.', 'Color', mygray(v1i,:), 'MarkerSize', 8);
end
%legend('V_1 = 1');
ylabel('Mean firing rate (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdirintersect, 12, [6, 2.5], 1);
