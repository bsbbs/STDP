% Testing the kernels
plotdir = 'C:\Users\Bo\NYU Langone Health Dropbox\Shen Bo\Bo Shen Working files\STDP_Project\SyncNMDA_1600';
Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
load(Ntwkfile);
load(fullfile(plotdir, 'RealtimeMonitor_Event3.mat'));
Ntwk.Synapse.tauNMDA.rise = 2; % ms, Feldmeyer et al., 2002
Ntwk.Synapse.tauNMDA.decay = 26; % Feldmeyer et al., 2002
dt = .1; % ms, time precision for simulation, in unit of second
dur = 500; % ms
% Time vector
timevec = [dt:dt:dur]';
tsteps = numel(timevec);
a = 0.062; % mV-1
b = 3.57; % mM
Mg2 = 1; % nM
Ntwk.Synapse.gbarE = 3.25; % original value .14 nS
Ntwk.Synapse.gbarI = 31.5; % original value .35 nS
ExctV = Ntwk.VL;
ExctV2 = Ntwk.VL;
InhbtV = Ntwk.VL;
ExctRefraction = 0;
ExctRefraction2 = 0;
InhbtRefraction = 0;
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
% - synaptic conductance (presynaptic-activity dependent)
ExctgAMPA = 0;
ExctgGABA = 0;
InhbtgAMPA = 0;
InhbtgI = 0;
xi = 0;
gNMDAi = 0;
gNMDA = 0;
Output = [];
for t = 1:tsteps
    % input spikes
    if t*dt == 5 || (mod(t*dt, 5) == 0 && t*dt >= 300 && t*dt < 400)
        InputSpikes = 1;
    else
        InputSpikes = 0;
    end
    % Synaptic plasticity
    WEI = 1;
    WIE = 1;
    WEE = 1;
    wInput = 1;
    % Synaptic activities
    ExctgAMPA = ExctgAMPA - ExctgAMPA/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarE*InputSpikes;
    ExctgGABA = ExctgGABA - ExctgGABA/Ntwk.Synapse.tauInhbt*dt ... % inhibitory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarI*InputSpikes;
    InhbtgAMPA = InhbtgAMPA - InhbtgAMPA/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Inhbt neurons
        + Ntwk.Synapse.gbarE*InputSpikes;
    xi = xi + (-xi/Ntwk.Synapse.tauNMDA.rise)*dt + InputSpikes;
    xi(xi < 0) = 0;
    gNMDAi = gNMDAi + (-gNMDAi/Ntwk.Synapse.tauNMDA.decay + (1-gNMDAi)*xi)*dt;
    gNMDAi(gNMDAi < 0) = 0;
    gNMDA = Ntwk.Synapse.gbarE*wInput*gNMDAi;
    
    % Membrane potential change for Exct neurons receiving excitatory input
    INMDA = gNMDA./(1+exp(-a*ExctV/b)*Mg2).*(Ntwk.VE - ExctV);
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + 0*ExctgAMPA.*(Ntwk.VE - ExctV) + INMDA)/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    ExctRefraction = ExctRefraction - 1;
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    ExctV(Espikes) = Ntwk.Vfire;
    ExctRefraction(Espikes) = refractionPeriod.E;

    % Membrane potential change for Exct neurons receiving inhibitory input
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV2) + ExctgGABA.*(Ntwk.VI - ExctV2))/Ntwk.Exct.Cm*dt;
    ExctV2 = ExctV2 + dV;
    ExctV2(ExctRefraction2>0) = Ntwk.Vreset;
    ExctRefraction2 = ExctRefraction2 - 1;
    % Exct neurons fire
    Espikes2 = ExctV2 > Ntwk.Vth;
    ExctV(Espikes2) = Ntwk.Vfire;
    ExctRefraction2(Espikes2) = refractionPeriod.E;
    
    % Membrane potential change for Inhbt neurons
    INMDA = gNMDA./(1+exp(-a*InhbtV/b)*Mg2).*(Ntwk.VE - InhbtV);
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + 0*InhbtgAMPA.*(Ntwk.VE - InhbtV) + INMDA)/Ntwk.Inhbt.Cm*dt;
    InhbtV = InhbtV + dV;
    InhbtV(InhbtRefraction>0) = Ntwk.Vreset;
    InhbtRefraction = InhbtRefraction - 1;
    % Inhbt neurons fire
    Ispikes = InhbtV > Ntwk.Vth;
    InhbtV(Ispikes) = Ntwk.Vfire;
    InhbtRefraction(Ispikes) = refractionPeriod.I;

    Output(t,:) = [ExctV, ExctV2, InhbtV]; 
end
%
h = figure; hold on;
filename = 'SynapseRestingAmplitudeNMDA';
lg = [];
lg(1) = plot(timevec, Output(:,1), 'k-', 'LineWidth',1);
lg(2) = plot(timevec, Output(:,2), 'k--', 'LineWidth',1);
lg(3) = plot(timevec, Output(:,3), 'r-', 'LineWidth',1);
legend(lg, {'Excitatory to E','Inhibitory to E','Excitatory to I'}, 'location', 'best');
ylabel('Membrane potential (mV)');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [4, 2], 1);

%% Connection probability
Distance = 0:200;
Ntwk.AxonRange.EE = 130; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and E
Ntwk.AxonRange.EI = 100; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and I
Ntwk.AxonRange.IE = 97; % um, standard deviation of Gaussian decay of connection probability over somatic distance of I and E
Ntwk.CnnctProb.EE = .1; % the maximum connection probabability from E to E
Ntwk.CnnctProb.EI = .2; % the maximum connection probabability from E to I
Ntwk.CnnctProb.IE = .3; % the maximum connection probabability from I to E


%% Single spike

%%
h = figure; hold on;
filename = 'LIF_demo_AP';
plot(timevec,Smpl.ExctV(:,1), 'k-');
ylim([-70, -45]);
xlim([1140, 1250]);
xticks([1150:50:1250]);
xticklabels({'0','50','100'});
ylabel('Membrane potential (mV)');
xlabel('Time (ms)');
plot(timevec,Smpl.InhbtV(:,1), 'r-');
mysavefig(h, filename, plotdir, 12, [2, 1.6], 1);
%% Input sequences
load(fullfile(plotdir, 'Seq.mat'));
leftt = evs(120,1);
rightt = evs(130,1);
Ntwk.Input.Source = 2;
h = figure;
filename = 'InputSequences';
subplot(3,1,1); hold on;
for ii = Ntwk.Input.Source:-1:1
    plot(time/1000, Seq(:,ii)*.9+(ii), '-', "Color", OKeeffe(ii,:), 'LineWidth', 1);
end
title('Input signal');
xlabel('Time (ms)');
ylabel('Channel');
axis([leftt-.5 rightt+.5 .5, Ntwk.Input.Source*1.45]); % Adjust the axis for better visualization
% ylim([.5, Ntwk.Input.Source*1.45]);
yticks([1:Ntwk.Input.Source]);
mysavefig(h, filename, plotdir, 12, [3,6], 1);
%% Dynamics of NMDA receptors
% rising and decay captured by two coupled linear differentail equations (Wilson and Bower, 1989; Destexhe et al. 1994)
Ntwk.Synapse.tauNMDA.rise = 2; % ms, Feldmeyer et al., 2002
Ntwk.Synapse.tauNMDA.decay = 26; % Feldmeyer et al., 2002
dt = .1; % ms, time precision for simulation, in unit of second
dur = 200; % ms
% Time vector
timevec = [dt:dt:dur]';
tsteps = numel(timevec);
wInput = 1;
gNMDAi = 0;
xi = 0; % the intermediate variable for the double exponential implementation
Output = [];
a = 0.062; % mV-1
b = 3.57; % mM
Mg2 = 1; % nM
V = -70:0;
Mgblock = 1./(1+exp(-a.*V./b).*Mg2);
h = figure; plot(V, Mgblock, '-');
ExctV = Ntwk.VL;
ExctRefraction = 0;
refractionPeriod.E = Ntwk.Exct.tauREF/dt;
refractionPeriod.I = Ntwk.Inhbt.tauREF/dt;
for t = 1:tsteps
    % input spikes
    if mod(t*dt, 5) == 0 && t*dt < 100 % t*dt == 5 || t*dt == 10 || t*dt == 15 || t*dt == 20
        InputSpikes = 1;
    else
        InputSpikes = 0;
    end
    dxi = (-xi/Ntwk.Synapse.tauNMDA.rise)*dt + InputSpikes;
    xi = xi + dxi;
    xi(xi < 0) = 0;
    gNMDAi = gNMDAi + (-gNMDAi/Ntwk.Synapse.tauNMDA.decay + (1-gNMDAi)*xi)*dt;
    gNMDAi(gNMDAi < 0) = 0;
    gNMDA = Ntwk.Synapse.gbarE*wInput*gNMDAi;
    INMDA = gNMDA./(1+exp(-a*ExctV/b)*Mg2).*(Ntwk.VE - ExctV);

    % Membrane potential change for Exct neurons receiving excitatory input
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + INMDA)/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    ExctRefraction = ExctRefraction - 1;
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    ExctV(Espikes) = Ntwk.Vfire;
    ExctRefraction(Espikes) = refractionPeriod.E;

    Output(t,:) = [INMDA, gNMDAi, xi, ExctV]; 
end
h = figure;
filename = 'NMDA_DoubleExpo';
subplot(1,3,1); hold on;
lg = [];
lg(1) = plot(timevec, Output(:,2), 'k--', 'LineWidth',1);
lg(2) = plot(timevec, Output(:,3), 'k:', 'LineWidth',1);
legend(lg, {'g_{NMDAi}','xi'});
ylabel('Activities');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [6, 2], 1);
subplot(1,3,2); hold on;
plot(timevec, Output(:,1), 'k-', 'LineWidth',1);
ylabel('Current (pA)');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [6, 2], 1);
subplot(1,3,3); hold on;
plot(timevec, Output(:,4), 'r-', 'LineWidth',1);
ylabel('Membrane potential');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [6, 2], 1);