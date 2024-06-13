% Testing the kernels
dt = .1; % ms, time precision for simulation, in unit of second
dur = 100; % ms
% Time vector
timevec = [dt:dt:dur]';
tsteps = numel(timevec);

Ntwk.Synapse.gbarE = 4.7; % original value .14 nS
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

    % Membrane potential change for Exct neurons receiving excitatory input
    dV = (Ntwk.Exct.gL*(Ntwk.VL - ExctV) + ExctgE.*(Ntwk.VE - ExctV))/Ntwk.Exct.Cm*dt;
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
    dV = (Ntwk.Inhbt.gL*(Ntwk.VL - InhbtV) + InhbtgE.*(Ntwk.VE - InhbtV))/Ntwk.Inhbt.Cm*dt;
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
ylabel('Time (ms)');
mysavefig(h, filename, plotdir, 12, [2, 2], 1);
