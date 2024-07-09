%% Evaluation of the computation of the circuit
%% Define the input and output directories
DefineIO;
Ntrial = 1600;
ProjectName = sprintf('SyncNMDA_%i', Ntrial);
plotdir = fullfile(Projdir, ProjectName);
evi = 0;
plotdirintersect = fullfile(plotdir, sprintf('Intersect%i', evi));
if ~exist(plotdirintersect,'dir')
    mkdir(plotdirintersect);
end
%% load Ntwk file and trained weights
Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
load(Ntwkfile);
filename = sprintf('RealtimeMonitor_Event%i', evi);
load(fullfile(plotdir, [filename, '.mat']));
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

duration = 40000*dt; % ms
time = [dt:dt:duration]';
timesteps = numel(time);
%% test the network
% Initializing the network status at t=0
% - membrane potentials
ExctV = Ntwk.VL*gpuArray.ones(Ntwk.Exct.N,1);
InhbtV = Ntwk.VL*gpuArray.ones(Ntwk.Inhbt.N,1);
ExctRefraction = gpuArray.zeros(Ntwk.Exct.N,1);
InhbtRefraction = gpuArray.zeros(Ntwk.Inhbt.N,1);
% - synaptic conductance (presynaptic-activity dependent)
ExctgAMPA = gpuArray.zeros(Ntwk.Exct.N, 1);
ExctgNMDA = gpuArray.zeros(Ntwk.Exct.N, 1);
xE = gpuArray.zeros(Ntwk.Exct.N, 1);
ExctgGABA = gpuArray.zeros(Ntwk.Exct.N, 1);
InhbtgAMPA = gpuArray.zeros(Ntwk.Inhbt.N, 1);
InhbtgNMDA = gpuArray.zeros(Ntwk.Inhbt.N, 1);
xI = gpuArray.zeros(Ntwk.Inhbt.N, 1);
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
bankwidth = max([Ntwk.Delay.EE, Ntwk.Delay.EE, Ntwk.Delay.IE])/dt;
SpkTrnE = nan(Ntwk.Exct.N+1,bankwidth); tickE = 0;
SpkTrnI = nan(Ntwk.Inhbt.N+1,bankwidth); tickI = 0;
% Output variables
MeanFR = [];
timevec = [];
for t = 1:(timesteps-1)
    % input spikes
    InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
    % InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(t,1);
    % InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(t,2);
    if t*dt <= 2000
        InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*1;
        InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*0;
    else
        InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*0;
        InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*0;
    end
    % Synaptic activities
    % AMPA on excitatory neurons
    ExctgAMPA = ExctgAMPA - ExctgAMPA/Ntwk.Synapse.tauExct*dt; % excitatory synaptic conductance on Exct neurons
    if any(InputSpikes)
        ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*Ntwk.Cnnct_Input.*Ntwk.wInput*InputSpikes;
    end
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EE/dt);
    if ~isempty(DlyEspikes)
        ExctgAMPA = ExctgAMPA + Ntwk.Synapse.gbarE*Ntwk.Cnnct_EE.*WEE*DlyEspikes;
    end
    % NMDA on excitatory neurons
    Delta = 0;
    if any(InputSpikes)
        Delta = Delta + Ntwk.Synapse.gbarE*Ntwk.Cnnct_Input.*Ntwk.wInput*InputSpikes;
    end
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EE/dt);
    if ~isempty(DlyEspikes)
        Delta = Delta + Ntwk.Synapse.gbarE*Ntwk.Cnnct_EE.*WEE*DlyEspikes;
    end
    xE = xE + (-xE/Ntwk.Synapse.tauNMDA.rise)*dt + Delta;
    ExctgNMDA = ExctgNMDA + (-ExctgNMDA/Ntwk.Synapse.tauNMDA.decay + xE)*dt;
    % GABA on excitatory neurons
    ExctgGABA = ExctgGABA - ExctgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Exct neurons
    DlyIspikes = SpkTrnI(2:end,SpkTrnI(1,:) == t - Ntwk.Delay.IE/dt);
    if ~isempty(DlyIspikes)
        ExctgGABA = ExctgGABA + Ntwk.Synapse.gbarI*Ntwk.Cnnct_IE.*WIE*DlyIspikes;
    end
    % AMPA on inhibitory neurons
    InhbtgAMPA = InhbtgAMPA - InhbtgAMPA/Ntwk.Synapse.tauExct*dt; % excitatory synaptic conductance on Inhbt neurons
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EI/dt);
    if ~isempty(DlyEspikes)
        InhbtgAMPA = InhbtgAMPA + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_EI.*WEI*DlyEspikes);
    end
    % NMDA on inhibitory neurons
    Delta = 0;
    DlyEspikes = SpkTrnE(2:end,SpkTrnE(1,:) == t - Ntwk.Delay.EI/dt);
    if ~isempty(DlyEspikes)
        Delta = Delta + Ntwk.Synapse.gbarE*(Ntwk.Cnnct_EI.*WEI*DlyEspikes);
    end
    xI = xI + (-xI/Ntwk.Synapse.tauNMDA.rise)*dt + Delta;
    InhbtgNMDA = InhbtgNMDA + (-InhbtgNMDA/Ntwk.Synapse.tauNMDA.decay + xI)*dt;
    % GABA on excitatory neurons
    ExctgGABA = ExctgGABA - ExctgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Exct neurons
    DlyIspikes = SpkTrnI(2:end,SpkTrnI(1,:) == t - Ntwk.Delay.IE/dt);
    if ~isempty(DlyIspikes)
        ExctgGABA = ExctgGABA + Ntwk.Synapse.gbarI*Ntwk.Cnnct_IE.*WIE*DlyIspikes;
    end
    % GABA on inhibitory neurons
    InhbtgGABA = InhbtgGABA - InhbtgGABA/Ntwk.Synapse.tauInhbt*dt; % inhibitory synaptic conductance on Inhbt neurons

    % Updating OU noise
    ExctNoise = ExctNoise + ((Ib - ExctNoise)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Exct.N,1)*Ntwk.Noise.sgm)*dt;
    InhbtNosie = InhbtNosie + ((Ib - InhbtNosie)/Ntwk.Noise.tauN + gpuArray.randn(Ntwk.Inhbt.N,1)*Ntwk.Noise.sgm)*dt;
    
    % Membrane potential change for Exct neurons
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + (ExctgAMPA + ExctgNMDA).*(Ntwk.VE - ExctV) + ExctgGABA.*(Ntwk.VI - ExctV) + ExctNoise)/Ntwk.Exct.Cm*dt;
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
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + (InhbtgAMPA + InhbtgNMDA).*(Ntwk.VE - InhbtV) + InhbtgGABA.*(Ntwk.VI - InhbtV) + InhbtNosie)/Ntwk.Inhbt.Cm*dt;
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
    MeanFR(t, 1) = sum(Espikes(TuneMask1))/sum(TuneMask1);
    MeanFR(t, 2) = sum(Espikes(TuneMask2))/sum(TuneMask2);
    timevec = [timevec; t*dt];

%     if any(t == smplonsets - 100/dt)
%         evi = find(t == smplonsets - 100/dt);
%         filename = sprintf('RealtimeMonitor_Event%i', evi);
%         % Prepare the video file
%         profiles = VideoWriter.getProfiles();
%         disp(profiles);
%         writerObj1 = VideoWriter(fullfile(plotdir, filename), 'Motion JPEG AVI');
%         writerObj1.FrameRate = 100; % Adjust frame rate as needed
%         writerObj1.Quality = 95;   % Set quality to maximum for best results (only for MPEG-4)
%         open(writerObj1);
%         timevec = [];
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
%         timevec = [timevec; t*dt];
%     end
    if mod(t*dt, 10) == 0
        fprintf('.');
        if mod(t*dt,100)==0
            fprintf('%1.4f s\n', t*dt/1000);
        end
    end
end
fprintf('Successfully completed\n');
%
h = figure; hold on;
filename = 'MeanFRDynmc';
[FR, timep] = PSTH(MeanFR, dt);
plot(timep, FR(:,1),'-', 'Color', OKeeffe(1,:), 'LineWidth', 2);
plot(timep, FR(:,2),'-', 'Color', OKeeffe(2,:), 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
mysavefig(h, filename, plotdirintersect, 12, [4,2], 1);