Ntwkfile = fullfile(gnrloutdir, 'Ntwk200.mat');
if ~exist(Ntwkfile, 'file')
    %% The structure of the network
    Ntwk.XScale = 150; % um, range of the x axis
    Ntwk.YScale = 0; % um, range of the y axis
    Ntwk.Exct.Location = []; % the physical location of excitatory cells
    Ntwk.Exct.N = 200; % number of the excitarory cells
    Ntwk.Exct.AxonRange = 130;
    Ntwk.Inhbt.Location = [];
    Ntwk.Inhbt.N = 50;
    Ntwk.Inhbt.AxonRange = 97; % um, the standard deviation of axon physical connection range for interneurons
    Ntwk.AxonRange.EE = 130; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and E
    Ntwk.AxonRange.EI = 100; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and I
    Ntwk.AxonRange.IE = 97; % um, standard deviation of Gaussian decay of connection probability over somatic distance of I and E
    Ntwk.CnnctProb.EE = .1; % the maximum connection probabability from E to E
    Ntwk.CnnctProb.EI = .2; % the maximum connection probabability from E to I
    Ntwk.CnnctProb.IE = .3; % the maximum connection probabability from I to E
    %% define LIF model parameters
    Ntwk.VL = -70; % mV, resting potential
    Ntwk.Vth = -50; % mV, threshold potential (Wang 2002; Vogels et al., 2011)
    Ntwk.Vfire = Ntwk.Vth; % mV, spiking peak potential
    Ntwk.Vreset = -60; % mV, reset potential (Vogels, 2011)
    Ntwk.VE = 0; % mV, excitatory synaptic potential
    Ntwk.VI = -80; % mV, inhibitory synaptic potential (Vogels et al., 2011; Burkitt et al., 2004)
    Ntwk.Exct.Cm = .5*1000; % nF (multiplied 1000 to match the scale of ms), membrane capacity for pyramidal neurons (Wang 2002)
    Ntwk.Inhbt.Cm = .2*1000; % nF (multiplied 1000 to match the scale of ms), membrane capacity for interneurons (Wang 2002)
    Ntwk.Exct.gL = 25; % nS, membrane leaky conductance for pyramidal neurons (Wang 2002)
    Ntwk.Inhbt.gL = 20; % nS, membrane leaky conductance for interneurons (Wang 2002)
    Ntwk.Exct.taum = Ntwk.Exct.Cm/Ntwk.Exct.gL; % 20 ms, membrane time constant of pyramidal neurons (Wang, 2002; Vogels et al., 2011; Burkitt et al., 2004)
    Ntwk.Inhbt.taum = Ntwk.Inhbt.Cm/Ntwk.Inhbt.gL; % 10 ms, membrane time constant of interneurons (Wang 2002)
    Ntwk.Exct.tauREF = 2; % ms, refractory period for pyramidal neurons (Wang 2002)
    Ntwk.Inhbt.tauREF = 1; % ms, refractory period for interneurons (Wang 2002)
    Ntwk.Delay.EE = 1.7; % ms, delay of presynaptic action potential to the PSP response onset (Campagnola et al., 2022am)
    Ntwk.Delay.EI = 1.3; % ms
    Ntwk.Delay.IE = 1; % ms
    Ntwk.Synapse.tauExct = 5; % ms, the decay time constant of EPSC, like AMPA (Destexhe et al., 1994; Gabbiani et al., 1994)
    Ntwk.Synapse.tauInhbt = 6; % ms, the decay time constant of IPSC, like GABAa (Destexhe & Paré, 1999)
    Ntwk.Synapse.NMDA.taurise = 2; % ms, Feldmeyer et al., 2002
    Ntwk.Synapse.NMDA.taudecay = 26; % Feldmeyer et al., 2002
    Ntwk.Synapse.NMDA.a = 0.062; % mV-1
    Ntwk.Synapse.NMDA.b = 3.57; % mM
    Ntwk.Synapse.NMDA.Mg2 = 1; % nM
    % rising time of those synapses were assumed as instant
    % Ntwk.Synapse.gbarE = .140; % nS, the maximum excitatory synaptic conductance, like AMPA (Vogels et al., 2011)
    % Ntwk.Synapse.gbarI = .350; % nS, the maximum inhibitory synaptic conductance, like GABA (Vogels et al., 2011)
    Ntwk.Synapse.gbarE = 3.25*4; %
    Ntwk.Synapse.gbarI = 31.5*4; %
    %% Spike-timing dependent plasticity
    % time constant of synaptic plasticity for pre-post and post-pre kernels
    Ntwk.A = .1; % the maximum amplitude change on the synapctic plasticity for each pair of spikes
    % changed from .01 for a faster equilibrium
    Ntwk.Synapse.tau_prepost = 20; % ms
    Ntwk.STDP.tau_postpre = 20; % ms
    Ntwk.ExctSTDP.tau_prepost = 20; % ms
    Ntwk.ExctSTDP.sign_prepost = 1;
    Ntwk.ExctSTDP.tau_postpre = 20; % ms
    Ntwk.ExctSTDP.sign_postpre = -1;
    Ntwk.InhbtSTDP.tau_prepost = 20; % ms
    Ntwk.InhbtSTDP.sign_prepost = 1;
    Ntwk.InhbtSTDP.tau_postpre = 20; % ms
    Ntwk.InhbtSTDP.sign_postpre = 1;
    Ntwk.ExctSTDP.eta = 1e-4; % the learning rate (updating speed) of synaptic connections
    Ntwk.ExctSTDP.intercept_pre = 0; % depression if no post spikes
    Ntwk.ExctSTDP.intercept_post = 0; % depression if no pre spikes
    Ntwk.InhbtSTDP.eta = 1e-4; % the learning rate (updating speed) of synaptic connections
    Ntwk.InhbtSTDP.rho = 5; % Hz, depression constant
    Ntwk.InhbtSTDP.intercept_pre = -2*Ntwk.InhbtSTDP.rho/1000*Ntwk.InhbtSTDP.tau_prepost; % depression if no post spikes
    Ntwk.InhbtSTDP.intercept_post = 0; % depression if no pre spikes
    %% Other parameters
    Ntwk.Noise.tauN = 2; % ms, time constant for the OU process of noise
    Ntwk.Noise.sgm = 30; % amplitude of noise

    %% Physical location of the cells, assuming located on the same layer, thus
    % spreading the 1-D space
    % locations of excitatory cells
    rng(2002);
    xE = -Ntwk.XScale + rand(Ntwk.Exct.N,1)*Ntwk.XScale*2;
    [xE, ~] = sort(xE);
    yE = Ntwk.YScale*ones(Ntwk.Exct.N,1)/2; % assume all neurons approximately on the same layer
    Ntwk.Exct.Location = [xE, yE];
    
    % locations of Inhibitory cells
    rng(2004);
    xI = -Ntwk.XScale + rand(Ntwk.Inhbt.N,1)*Ntwk.XScale*2;
    [xI, ~] = sort(xI);
    yI = -Ntwk.YScale*ones(Ntwk.Inhbt.N,1)/2;
    Ntwk.Inhbt.Location = [xI, yI];
    clear xE yE xI yI;

    %% Establish possible synaptic connections, independent from synaptic weights
    % possible synaptic connection from E -> I
    [XE, XI] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1)); % rows represent Inhbt and columns represent Exct
    [YE, YI] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Inhbt.Location(:,2));
    DstcEI = sqrt((XE - XI).^2 + (YE - YI).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
    % or periodic
    % DstcEI = sqrt(min(Ntwk.XScale*2-abs(XE - XI), abs(XE - XI)).^2 + min(Ntwk.YScale*2 - abs(YE - YI), abs(YE - YI)).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
    p_EI = Ntwk.CnnctProb.EI*exp(-.5*(DstcEI/Ntwk.AxonRange.EI).^2); %1 - (1 - exp(-DstcEI/Ntwk.Exct.AxonRange)).^4; % probability of physical connection from E to I based on distance
    rng(2007);
    Ntwk.Cnnct_EI = p_EI >= rand(size(p_EI)); % connections from E to I, 0 or 1
    clear XE XI YE YI DstcEI p_EI;
    % possible synaptic connection from I -> E
    [XI, XE] = meshgrid(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent Exct and columns represent Inhbt
    [YI, YE] = meshgrid(Ntwk.Inhbt.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcIE = sqrt((XI - XE).^2 + (YI - YE).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
    % or periodic
    % DstcIE = sqrt(min(Ntwk.XScale*2 - abs(XI - XE), abs(XI - XE)).^2 + min(Ntwk.YScale*2 - abs(YI - YE), abs(YI - YE)).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
    p_IE = Ntwk.CnnctProb.IE*exp(-.5*(DstcIE/Ntwk.AxonRange.IE).^2); % probability of connection from I to E based on distance
    rng(2008);
    Ntwk.Cnnct_IE = p_IE >= rand(size(p_IE)); % connections from I to E, 0 or 1
    clear XE XI YE YI DstcIE p_IE;
    % possible synaptic connection from E -> E
    [XE1, XE2] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1));
    [YE1, YE2] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcEE = sqrt((XE1 - XE2).^2 + (YE1 - YE2).^2); %  Euclidean distance between each pair of Exct neurons
    % or periodic
    % DstcEE = sqrt(min(Ntwk.XScale*2-abs(XE1 - XE2), abs(XE1 - XE2)).^2 + min(Ntwk.YScale*2 - abs(YE1 - YE2), abs(YE1 - YE2)).^2); %  Euclidean distance between each pair of Exct neurons
    p_EE = Ntwk.CnnctProb.EE*exp(-.5*(DstcEE/Ntwk.AxonRange.EE).^2); % probability of connection from E to E based on distance
    rng(2009);
    Ntwk.Cnnct_EE = p_EE >= rand(size(p_EE)); % connections from E to E, 0 or 1
    clear XE1 XE2 YE1 YE2 DstcEE p_EE;
    % possible synaptic connection from I -> I was not assumed since we do
    % not assume disinhibition in the current interneuronal type (PV like neurons)
    %% Visualization
    downsmpl = 5;
    Ntwk.Smpl.E = [66, 134];
    Ntwk.Smpl.I = [17, 33];
    if visualize
        h = figure;
        filename = 'Network structure';
        hold on;
        plot(Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,1), Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,2), 'k^', 'MarkerSize', 4, 'MarkerFaceColor','auto');
        plot(Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,1),...
            Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,2)-1, '.', 'Color', 'r', 'MarkerSize', 10);
        xlabel('x (\mum)');
        ylabel(' ');
        % axis equal;
        xlim([-Ntwk.XScale, Ntwk.XScale]);
        ylim([-1.5, .5]);
        mysavefig(h, filename, gnrloutdir, 14, [4, 2], 2);

        filename = 'Network structure - Example EIs';
        for ei = Ntwk.Smpl.E
            EIs = find(Ntwk.Cnnct_EI(:, ei));
            for i = EIs'
                lgd1 = plot([Ntwk.Exct.Location(ei,1), Ntwk.Inhbt.Location(i,1)],...
                    [Ntwk.Exct.Location(ei,2), Ntwk.Inhbt.Location(i,2)-1],...
                    '-', 'Color', 'b', 'LineWidth',2);
            end
        end
        %plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', Mycolors(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);

        for ii = Ntwk.Smpl.I
            IEs = find(Ntwk.Cnnct_IE(:, ii));
            for i = IEs'
                lgd2 = plot([Ntwk.Exct.Location(i,1), Ntwk.Inhbt.Location(ii,1)],...
                    [Ntwk.Exct.Location(i,2), Ntwk.Inhbt.Location(ii,2)-1],...
                    '--', 'Color', 'r', 'LineWidth',2);
            end
        end
        % plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', Mycolors(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        legend([lgd1, lgd2], {'E to I','I to E'});
        mysavefig(h, filename, gnrloutdir, 14, [4, 2], 2);
        clear EIs IEs i;
    end
    %% Synaptic weights for those possible connections defined above
    % initial weights
    rng(2010);
    Ntwk.wEI_initial = .08*rand(size(Ntwk.Cnnct_EI)).*Ntwk.Cnnct_EI; % synaptic weight from E to I, weak initial connections from E to I
    Ntwk.wIE_initial = .08*rand(size(Ntwk.Cnnct_IE)).*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST
    Ntwk.wEE_initial = .08*rand(size(Ntwk.Cnnct_EE)).*Ntwk.Cnnct_EE; % synaptic weight from E to E, weak initial connections of self-excitation
    %% visualization
    if visualize
        h = figure;
        filename = 'wEE_initial';
        imagesc(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1), Ntwk.wEE_initial);
        caxis([0, .1]);
        colormap(bluewhitered(256));
        c = colorbar;
        c.Label.String = 'Initial synaptic weights';
        c.Location = 'northoutside';
        xlabel('X-axis of E (\mum)');
        ylabel('X-axis of E (\mum)');
        axis equal;
        xlim([-Ntwk.XScale, Ntwk.XScale]);
        ylim([-Ntwk.XScale, Ntwk.XScale]);
        mysavefig(h, filename, gnrloutdir, 12, [4, 2]); % [4, 3] % [3.2, 3.6]

        h = figure;
        filename = 'wEI_initial';
        imagesc(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1), Ntwk.wEI_initial);
        caxis([0, .1]);
        colormap(bluewhitered(256));
        c = colorbar;
        c.Label.String = 'Initial synaptic weights';
        c.Location = 'northoutside';
        xlabel('X-axis of E (\mum)');
        ylabel('X-axis of I (\mum)');
        axis equal;
        xlim([-Ntwk.XScale, Ntwk.XScale]);
        ylim([-Ntwk.XScale, Ntwk.XScale]);
        mysavefig(h, filename, gnrloutdir, 12, [4, 2]);

        h = figure;
        filename = 'wIE_initial';
        imagesc(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1), Ntwk.wIE_initial);
        caxis([0, .1]);
        colormap(bluewhitered(256));
        c = colorbar;
        c.Label.String = 'Initial synaptic weights';
        c.Location = 'northoutside';
        ylabel('X-axis of E (\mum)');
        xlabel('X-axis of I (\mum)');
        axis equal;
        xlim([-Ntwk.XScale, Ntwk.XScale]);
        ylim([-Ntwk.XScale, Ntwk.XScale]);
        mysavefig(h, filename, gnrloutdir, 12, [4, 2]);
    end

    %% Inputs of two sources, assume inputs only intervene Exct neurons
    %% Define input connection matrix
    Ntwk.Input.Source = 2; % Number of input source(s)
    Ntwk.Input.N = 100*Ntwk.Input.Source; % number of the input (excitarory) neurons
    Ntwk.Input.XLoc = [-50, 50]; % The x-axis locations of the two input sources.
    Ntwk.AxonRange.Input = 50; % 130; % um, for each input neuron, the axon connection to the neighbouring E neuron is assumed the same as Gaussian distance delay as E to E
    Ntwk.CnnctProb.Input = .2; % The connection probability from input neurons are set as the same as E to E
    Ntwk.Input.on = .5; % input intervene 50 % of the neurons nearby
    Ntwk.Input.spikeRate = 100/1000; % firing rate per second
    tmp = repmat(Ntwk.Input.XLoc, Ntwk.Input.N/Ntwk.Input.Source, 1);
    xE = tmp(:); 
    yE = ones(size(xE))*Ntwk.YScale; 
    tmp = repmat([1, 2], Ntwk.Input.N/Ntwk.Input.Source, 1);
    Origins = tmp(:);
    Ntwk.Input.Location = [xE, yE];
    Ntwk.Input.Origins = Origins;
    clear xE yE Origins;
    %% connections to the local excitatory neurons
    [XInput, XE] = meshgrid(Ntwk.Input.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent local and columns represent Input projection
    [YInput, YE] = meshgrid(Ntwk.Input.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcInput = sqrt((XInput - XE).^2 + (YInput - YE).^2); % Euclidean distance between each pair of neurons
    % or periodic
    % DstcInput = sqrt(min(Ntwk.XScale*2 - abs(XInput - XE), abs(XInput - XE)).^2 + min(Ntwk.YScale*2 - abs(YInput - YE), abs(YInput - YE)).^2); % Euclidean distance between each pair of neurons
    p_Input = Ntwk.CnnctProb.Input*exp(-.5*(DstcInput/Ntwk.AxonRange.Input).^2); % probability of physical connection based on distance
    rng(2097);
    Ntwk.Cnnct_Input = p_Input >= rand(size(p_Input));
    for i = 1:Ntwk.Input.Source % set on and off neurons
        rng(2036+i);
        tmp = Ntwk.Cnnct_Input(:,Ntwk.Input.Origins == i);
        Ntwk.Cnnct_Input(:,Ntwk.Input.Origins == i) = repmat((Ntwk.Input.on >= rand(Ntwk.Exct.N,1)), 1, Ntwk.Input.N/Ntwk.Input.Source) .* tmp; % projections from input, 0 or 1
    end
    Ntwk.wInput = rand(size(Ntwk.Cnnct_Input)).*Ntwk.Cnnct_Input;
    clear XInput YInput XE YE DstcInput p_Input;
    save(Ntwkfile, 'Ntwk');
    %%
    if visualize
        h = figure;
        filename = 'InputTuning';
        hold on;
        plot(Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,1), Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,2), 'k^', 'MarkerSize', 4, 'MarkerFaceColor','auto');
        plot(Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,1),...
            Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,2)-1, '.', 'Color', 'r', 'MarkerSize', 10);
        for Surci = 1:Ntwk.Input.Source
            plot(Ntwk.Input.XLoc(Surci), 1, '.', 'Color', OKeeffe(Surci,:), 'MarkerSize', 10);
        end
        % mapping connections
        for ei = Ntwk.Smpl.E % from E to I
            EIs = find(Ntwk.Cnnct_EI(:, ei));
            for i = EIs'
                lgd1 = plot([Ntwk.Exct.Location(ei,1), Ntwk.Inhbt.Location(i,1)],...
                    [Ntwk.Exct.Location(ei,2), Ntwk.Inhbt.Location(i,2)-1],...
                    '-', 'Color', 'b', 'LineWidth',2);
            end
        end
        for ii = Ntwk.Smpl.I % from I to E
            IEs = find(Ntwk.Cnnct_IE(:, ii));
            for i = IEs'
                lgd2 = plot([Ntwk.Exct.Location(i,1), Ntwk.Inhbt.Location(ii,1)],...
                    [Ntwk.Exct.Location(i,2), Ntwk.Inhbt.Location(ii,2)-1],...
                    '--', 'Color', 'r', 'LineWidth',2);
            end
        end
        Ntwk.Smpl.Input = [70, 130];
        Surci = 0;
        for Inpti = Ntwk.Smpl.Input % from Inputs to E
            Surci = Surci + 1;
            InptEs = find(Ntwk.Cnnct_Input(:, Inpti));
            for i = InptEs'
                lgd2 = plot([Ntwk.Input.Location(Inpti,1), Ntwk.Exct.Location(i,1)],...
                    [Ntwk.Input.Location(Inpti,2)+1, Ntwk.Exct.Location(i,2)],...
                    '-', 'Color', OKeeffe(Surci,:), 'LineWidth',2);
            end
        end
        xlabel('x (\mum)');
        xlim([-Ntwk.XScale, Ntwk.XScale]);
        ylim([-1.5, 1.2]);
        mysavefig(h, filename, gnrloutdir, 14, [4, 2], 2);
    end
else
    load(Ntwkfile, 'Ntwk');
end
