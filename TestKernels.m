% Testing the kernels
dt = .1; % ms, time precision for simulation, in unit of second
dur = 100; % ms
% Time vector
timevec = [dt:dt:dur]';
tsteps = numel(timevec);

Ntwk.Synapse.gbarE = 3.75; % original value .14 nS
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
ExctgE = 0;
ExctgI = 0;
InhbtgE = 0;
InhbtgI = 0;
x = 0;
gNMDA = 0;
Output = [];
for t = 1:tsteps
    % input spikes
    if t*dt == 5
        InputSpikes = 1;
    else
        InputSpikes = 0;
    end
    % Synaptic plasticity
    WEI = 1;
    WIE = 1;
    WEE = 1;
    % Synaptic activities
    ExctgE = ExctgE - ExctgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarE*InputSpikes;
    ExctgI = ExctgI - ExctgI/Ntwk.Synapse.tauInhbt*dt ... % inhibitory synaptic conductance on Exct neurons
        + Ntwk.Synapse.gbarI*InputSpikes;
    InhbtgE = InhbtgE - InhbtgE/Ntwk.Synapse.tauExct*dt ... % excitatory synaptic conductance on Inhbt neurons
        + Ntwk.Synapse.gbarE*InputSpikes;
    dx = (-x/Ntwk.Synapse.tauNMDA.rise)*dt + Ntwk.Synapse.gbarE/18*InputSpikes;
    x = x + dx;
    x(x < 0) = 0;
    dgNMDA = (-gNMDA/Ntwk.Synapse.tauNMDA.decay + x)*dt;
    gNMDA = gNMDA + dgNMDA;
    gNMDA(gNMDA < 0) = 0;
    
    % Membrane potential change for Exct neurons receiving excitatory input
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + (ExctgE + gNMDA).*(Ntwk.VE - ExctV))/Ntwk.Exct.Cm*dt;
    ExctV = ExctV + dV;
    ExctV(ExctRefraction>0) = Ntwk.Vreset;
    ExctRefraction = ExctRefraction - 1;
    % Exct neurons fire
    Espikes = ExctV > Ntwk.Vth;
    ExctV(Espikes) = Ntwk.Vfire;
    ExctRefraction(Espikes) = refractionPeriod.E;

    % Membrane potential change for Exct neurons receiving inhibitory input
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV2) + ExctgI.*(Ntwk.VI - ExctV2))/Ntwk.Exct.Cm*dt;
    ExctV2 = ExctV2 + dV;
    ExctV2(ExctRefraction2>0) = Ntwk.Vreset;
    ExctRefraction2 = ExctRefraction2 - 1;
    % Exct neurons fire
    Espikes2 = ExctV2 > Ntwk.Vth;
    ExctV(Espikes2) = Ntwk.Vfire;
    ExctRefraction2(Espikes2) = refractionPeriod.E;
    
    % Membrane potential change for Inhbt neurons
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + (InhbtgE + gNMDA).*(Ntwk.VE - InhbtV))/Ntwk.Inhbt.Cm*dt;
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
filename = 'SynapseRestingAmplitude';
lg = [];
lg(1) = plot(timevec, Output(:,1), 'k-', 'LineWidth',1);
lg(2) = plot(timevec, Output(:,2), 'k--', 'LineWidth',1);
lg(3) = plot(timevec, Output(:,3), 'r-', 'LineWidth',1);
legend(lg, {'Excitatory to E','Inhibitory to E','Excitatory to I'});
ylabel('Membrane potential');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [2, 2], 1);

%% Connection probability
Distance = 0:200;
Ntwk.AxonRange.EE = 130; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and E
Ntwk.AxonRange.EI = 100; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and I
Ntwk.AxonRange.IE = 97; % um, standard deviation of Gaussian decay of connection probability over somatic distance of I and E
Ntwk.CnnctProb.EE = .1; % the maximum connection probabability from E to E
Ntwk.CnnctProb.EI = .2; % the maximum connection probabability from E to I
Ntwk.CnnctProb.IE = .3; % the maximum connection probabability from I to E


%% Single spike
plotdir = 'C:\Users\Bo\Documents\STDP_Project\Sync_1600';
load(fullfile(plotdir, 'RealtimeMonitor_Event1.mat'));
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
Ntwk.Synapse.tauNMDA.rise = 3; % ms
Ntwk.Synapse.tauNMDA.decay = 100; % ms
dt = .1; % ms, time precision for simulation, in unit of second
dur = 100; % ms
% Time vector
timevec = [dt:dt:dur]';
tsteps = numel(timevec);
gNMDA = 0;
x = 0; % the intermediate variable for the double exponential implementation
Output = [];
for t = 1:tsteps
    % input spikes
    if t*dt == 5
        InputSpikes = 1;
    else
        InputSpikes = 0;
    end
    dx = (-x/Ntwk.Synapse.tauNMDA.rise)*dt + Ntwk.Synapse.gbarE*InputSpikes;
    x = x + dx;
    x(x < 0) = 0;
    dgNMDA = (-gNMDA/Ntwk.Synapse.tauNMDA.decay + x)*dt;
    gNMDA = gNMDA + dgNMDA;
    gNMDA(gNMDA < 0) = 0;
    Output(t,:) = [gNMDA, x]; 
end
h = figure; hold on;
filename = 'NMDA_DoubleExpo';
lg = [];
lg(1) = plot(timevec, Output(:,1), 'k-', 'LineWidth',1);
lg(2) = plot(timevec, Output(:,2), 'k--', 'LineWidth',1);
legend(lg, {'gNMDA','x'});
ylabel('Activities');
xlabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [2, 2], 1);
