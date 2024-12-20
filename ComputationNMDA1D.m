%% Evaluation of the computation of the circuit
%% Define the input and output directories
% DefineIO1D;
% dt = .1;
% Ntrial = 1600;
% ProjectName = sprintf('SyncNMDA1D100_%i', Ntrial);
% plotdir = fullfile(Projdir, ProjectName);
%% 
evi = 18; % intersectional check based on training duration from 0 to 18 period.
plotdirintersect = fullfile(subplotdir, sprintf('Intersect%i', evi));
if ~exist(plotdirintersect,'dir')
    mkdir(plotdirintersect);
    plotdirintersect_sim = fullfile(plotdirintersect, "SavedSimulations");
    mkdir(plotdirintersect_sim);
end
%% load Ntwk file and trained weights
Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
load(Ntwkfile);
Rsltfile = fullfile(subplotdir,'Rslts.mat');
load(Rsltfile);
if evi >=1
    filename = sprintf('RealtimeMonitor_Event%i', evi);
    load(fullfile(subplotdir, [filename, '.mat']));
else
    WEI = Ntwk.wEI_initial;
    WIE = Ntwk.wIE_initial;
    WEE = Ntwk.wEE_initial;
end
% EvalTuning1D(Ntwk,WEE,WEI,WIE,OKeeffe,plotdirintersect);
%% target the tuning neurons
Cnnct1 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 1), 2);
Cnnct2 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 2), 2);
Cnnct = sum(Ntwk.Cnnct_Input, 2);
InhbtCnnct = Ntwk.Cnnct_EI*Cnnct;
TuneMask1 = Cnnct1 & ~Cnnct2;
TuneMask2 = ~Cnnct1 & Cnnct2;
TuneMask = Cnnct > 3;
TuneMaskInhbt = InhbtCnnct > 200;
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

duration = 600; % ms
time = [dt:dt:duration]';
timesteps = numel(time);
% test the network
V1vec = 0:.5:2;
V2vec = flip(1./linspace(1/10,1,30)-1);
Modulation = nan(numel(V1vec),numel(V2vec), 2);
for v1i = 1:numel(V1vec)
    for v2i = 1:numel(V2vec)
        V1 = V1vec(v1i);
        V2 = V2vec(v2i);
        filename = sprintf('RealtimeMonitor_Event%i_testV1%1.1f_V2%1.1f', evi, V1, V2);
        if ~exist(fullfile(plotdirintersect_sim, [filename, '.mat']), 'file')
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
            ylabel('x (\mum)');
            ylim([-Ntwk.XScale, Ntwk.XScale]);
            % ylim([1, Ntwk.Exct.N]);
            % title('Raster plot of Exct/Inhbt neurons');
            ExctVec = 1:Ntwk.Exct.N;
            InhbtVec = [1:Ntwk.Inhbt.N]*4;
            % Output variables
            MeanFR = [];
            MeanFREI = [];
            timevec = [];
            for t = 1:(timesteps-1)
                % input spikes
                InputSpikes = gpuArray.rand(Ntwk.Input.N,1);
                % InputSpikes(Ntwk.Input.Origins == 1) = InputSpikes(Ntwk.Input.Origins == 1) < spikeProbability*Seq(t,1);
                % InputSpikes(Ntwk.Input.Origins == 2) = InputSpikes(Ntwk.Input.Origins == 2) < spikeProbability*Seq(t,2);
                if t*dt <= 550 && t*dt >= 50
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
                
                plot(time(t)*ones(sum(Espikes),1), Ntwk.Exct.Location(Espikes,1), 'k.', 'MarkerSize', 2);
                plot(time(t)*ones(sum(Ispikes),1), Ntwk.Inhbt.Location(Ispikes,1), 'r.', 'MarkerSize', 3);
                % plot(time(t)*ones(sum(Ispikes),1), InhbtVec(Ispikes), 'r.', 'MarkerSize', 3);
                % Example traces
                MeanFR(t, 1) = sum(Espikes(TuneMask1))/sum(TuneMask1);
                MeanFR(t, 2) = sum(Espikes(TuneMask2))/sum(TuneMask2);
                MeanFREI(t, 1) = sum(Espikes(TuneMask))/sum(TuneMask);
                MeanFREI(t, 2) = sum(Ispikes(TuneMaskInhbt))/sum(TuneMaskInhbt);
                timevec = [timevec; t*dt];
                if mod(t*dt, 10) == 0
                    fprintf('.');
                    if mod(t*dt,100)==0
                        fprintf('%1.4f s\n', t*dt/1000);
                    end
                end
            end
            fprintf('Successfully completed\n');
            mysavefig(h1, filename, plotdirintersect_sim, 12, [4, 4], 1);
            savefig(h1, fullfile(plotdirintersect_sim,[filename, '.fig']));
            save(fullfile(plotdirintersect_sim, [filename, '.mat']), 'timevec','MeanFR','MeanFREI');
        else
            load(fullfile(plotdirintersect_sim, [filename, '.mat']));
        end
        %
        h = figure; hold on;
        filename = sprintf('MeanFRDynmc_V1%1.0fV2%.1f', V1, V2);
        [FR, timep] = PSTH(MeanFR, dt);
        plot(timep, FR(:,1),'-', 'Color', OKeeffe(1,:), 'LineWidth', 2);
        plot(timep, FR(:,2),'-', 'Color', OKeeffe(2,:), 'LineWidth', 2);
        xlabel('Time (ms)');
        ylabel('Firing rate (Hz)');
        mysavefig(h, filename, plotdirintersect_sim, 12, [4, 2], 1);
        
        Modulation(v1i, v2i, :) = mean(FR(timep>200 & timep<1000,:),1);
    end
end
%%
h = figure; 
filename = 'ModulationEffect';
mygray = flip(gray(numel(V1vec)+1));
subplot(1,2,1);hold on; % direct effect
for v1i = 1:numel(V1vec)
    plot(V2vec, Modulation(v1i,:,2), 'k.', 'Color', mygray(v1i+1,:), 'MarkerSize', 8);
end
legend({'V_1 = 0', 'V_1 = 0.5','V_1 = 1.0','V_1 = 1.5','V_1 = 2.0'},'Location','best');
ylabel('Mean firing rate (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdirintersect, 12, [6, 2.5], 1);
subplot(1,2,2);hold on; % context effect
for v1i = 1:numel(V1vec)
    plot(V2vec, Modulation(v1i,:,1), 'k.', 'Color', mygray(v1i+1,:), 'MarkerSize', 8);
end
%legend('V_1 = 1');
ylabel('Mean firing rate (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdirintersect, 12, [6, 2.5], 1);
%% Model fit to the modulation effect
[V2mtx, V1mtx] = meshgrid(V2vec,1);
V2vecp = [0:.05:2];
[V2mtxp, V1mtxp] = meshgrid(V2vecp,1);
% [V2mtx, V1mtx] = meshgrid(V2vec,V1vec);
data = Modulation(2,:,1);
Comparisons = [];
Rslts = table('Size', [0 10], 'VariableTypes', {'uint8', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'},...
    'VariableNames', {'modeli', 'name', 'Mp', 'w1p', 'w2p', 'a', 'b1', 'b2', 'RSS', 'rsqrd'});
for modeli = 1:3
    if modeli == 1
        % divisive normalization
        modelname = 'DNM';
        RSSfun = @(params) DNM(params, data, V1mtx, V2mtx);
        % M, w1, w2
        PLB = [0, 0, 0]; % soft boundaries
        PUB = [1, 1, 1];
        x0 = rand(1, numel(PLB)).*(PUB - PLB) + PLB;
    elseif modeli == 2
        % linear subtraction
        modelname = 'LnSbtrct';
        RSSfun = @(params) LnSbtrct(params, data, V1mtx, V2mtx);
        % a, b1, b2
        PLB = [-5, 0, -5]; % soft boundaries
        PUB = [15, 15, 1];
        x0 = rand(1, numel(PLB)).*(PUB - PLB) + PLB;
    elseif modeli == 3
        % mixed
        modelname = 'Mixed';
        RSSfun = @(params) Mixed(params, data, V1mtx, V2mtx);
        % M, w1, w2, b
        PLB = [0, 0, 0, -10]; % soft boundaries
        PUB = [1, 1, 1, 10];
        x0 = rand(1, numel(PLB)).*(PUB - PLB) + PLB;
    end
    [fvalbest, ~, ~] = RSSfun(x0);
    fprintf('test succeeded\n');
    %%
    x0 = rand(1, numel(PLB)).*(PUB - PLB) + PLB;
    x = fmincon(RSSfun, x0);
    [RSS, rsqrd, Predicted] = eval([modelname, '(x, data, V1mtx, V2mtx)']);
    if modeli == 1
        new_row = table(modeli, {modelname}, x(1), x(2), x(3), nan,nan,nan, RSS, rsqrd, 'VariableNames', Rslts.Properties.VariableNames);
        tmp = (V1mtxp.^2)./(x(1) + x(2)*(V1mtxp.^2) + x(3)*(V2mtxp.^2));
        tmp(tmp<0) = 0;
    elseif modeli == 2
        new_row = table(modeli, {modelname}, nan,nan,nan, x(1), x(2), x(3), RSS, rsqrd, 'VariableNames', Rslts.Properties.VariableNames);
        tmp = x(1) + x(2)*V1mtxp + x(3)*V2mtxp;
        tmp(tmp<0) = 0;
    elseif modeli == 3
        new_row = table(modeli, {modelname}, x(1), x(2), x(3), nan,nan, x(4), RSS, rsqrd, 'VariableNames', Rslts.Properties.VariableNames);
        alpha = 1;
        tmp = (V1mtxp.^alpha)./(x(1) + x(2)*(V1mtxp.^alpha) + x(3)*(V2mtxp.^alpha)) + x(4)*(V2mtxp.^alpha);
        tmp(tmp<0) = 0;
    end
    Comparisons(modeli,:) = tmp;
    Rslts = [Rslts; new_row];
    writetable(Rslts, fullfile(plotdirintersect, 'ModelFitRslts.txt'), 'Delimiter', '\t');
    %%
    h = figure; hold on; % direct effect
    filename = sprintf('ModelPrediction_%s', modelname);
    mygray = flip(gray(numel(V1vec)+1));
    %for v1i = 2%1:numel(V1vec)
        plot(V2vec, data, 'k.', 'MarkerSize', 8);
        plot(V2vec, Predicted, 'k-', 'LineWidth', 1)
    %end
    legend(sprintf("r^2 = %0.3f", rsqrd), 'Location', 'best');
    % text(sprintf("r^2 = %0.3f", rsqrd));
    ylabel('Mean firing rate (Hz)');
    xlabel('V_2 (a.u.)');
    mysavefig(h, filename, plotdirintersect, 12, [3, 2.5], 1);
end


%% Model comparison
h = figure; hold on; % direct effect
filename = sprintf('ModelPrediction_Compare');
lg = [];
lg(1) = plot(V2vec, data, 'k.', 'MarkerSize', 8);
lg(2) = plot(V2vecp, Comparisons(1,:), '-', 'LineWidth', 1);
lg(3) = plot(V2vecp, Comparisons(2,:), '-', 'LineWidth', 1);
legend(lg, {'Data','Divisive normalization','Linear Subtraction'}, 'Location', 'best');
ylim([0,14]);
ylabel('Mean firing rate (Hz)');
xlabel('V_2 (a.u.)');
mysavefig(h, filename, plotdirintersect, 12, [3, 2.5], 1);

%%
function [RSS, rsqrd, Predicted] = DNM(params, data, V1mtx, V2mtx) 
M = params(1);
w1 = params(2);
w2 = params(3);
alpha = 2; % params(4);
Predicted = (V1mtx.^alpha)./(M + w1*(V1mtx.^alpha) + w2*(V2mtx.^alpha));
Predicted(Predicted<0) = 0;
% residual sum of squares
RSS = sum((data - Predicted).^2, 'all');
% r-squared
rsqrd = 1 - (RSS/sum((data - mean(data(:))).^2, 'all'));
end

function [RSS, rsqrd, Predicted] = LnSbtrct(params, data, V1mtx, V2mtx) 
a = params(1);
b1 = params(2);
b2 = params(3);
alpha = 1; % params(4);
Predicted = a + b1*V1mtx.^alpha + b2*V2mtx.^alpha;
Predicted(Predicted<0) = 0;
% residual sum of squares
RSS = sum((data - Predicted).^2, 'all');
% r-squared
rsqrd = 1 - (RSS/sum((data - mean(data(:))).^2, 'all'));
end

function [RSS, rsqrd, Predicted] = Mixed(params, data, V1mtx, V2mtx) 
M = params(1);
w1 = params(2);
w2 = params(3);
b = params(4);
alpha = 1; %params(4);
Predicted = (V1mtx.^alpha)./(M + w1*(V1mtx.^alpha) + w2*(V2mtx.^alpha)) + b*(V2mtx.^alpha);
Predicted(Predicted<0) = 0;
% residual sum of squares
RSS = sum((data - Predicted).^2, 'all');
% r-squared
rsqrd = 1 - (RSS/sum((data - mean(data(:))).^2, 'all'));
end
