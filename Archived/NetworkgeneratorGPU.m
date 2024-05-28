Ntwkfile = fullfile(gnrloutdir, 'NtkwGPU.mat');
if ~exist(Ntwkfile, 'file')
    %% The structure of the network
    Ntwk.Scale = 1000; % um, scale of the micro structure to test, assumed as one dimentional 
    Ntwk.Exct.Location = []; % the physical location of excitatory cells
    Ntwk.Exct.N = 8000; % number of the excitarory cells
    Ntwk.Exct.AxonRange = 150; % um, the standard deviation of axon physical connection range for pyramidal neurons
    Ntwk.Inhbt.Location = [];
    Ntwk.Inhbt.N = 2000;
    Ntwk.Inhbt.AxonRange = 50; % um, the standard deviation of axon physical connection range for interneurons
    
    %% Plasticity kernel
    % time constant of synaptic plasticity for pre-post and post-pre kernels
    Ntwk.A = .001; % the maximum amplitude change on the synapctic plasticity for each pair of spikes
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
    Ntwk.InhbtSTDP.intercept_pre = -2* Ntwk.InhbtSTDP.rho/1000* Ntwk.InhbtSTDP.tau_prepost; % depression if no post spikes
    Ntwk.InhbtSTDP.intercept_post = 0; % depression if no pre spikes
    %% define LIF model parameters
    Ntwk.VL = -60; % mV, resting potential
    Ntwk.Vth = -50; % mV, threshold potential (Wang 2002; Vogels et al., 2011)
    Ntwk.Vfire = 10; % mV, spiking peak potential
    Ntwk.Vreset = -55; % mV, reset potential (Wang 2002)
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
    Ntwk.Synapse.tauExct = 5; % 2; % ms, the decay time constant of EPSC, e.g., AMPAa (Wang 2002)
    Ntwk.Synapse.tauInhbt = 10; % 5; % ms, the decay time constant of IPSC, e.g., GABA (Wang 2002)
    Ntwk.Synapse.gbarE = .140; % nS, the maximum excitatory synaptic conductance, e.g., AMPAa (Vogels et al., 2011)
    Ntwk.Synapse.gbarI = .350; % nS, the maximum inhibitory synaptic conductance, e.g., GABA (Vogels et al., 2011)

    %% Other parameters
    Ntwk.Noise.tauN = 2; % ms, time constant for the OU process of noise
    Ntwk.Noise.sgm = 10; % amplitude of noise

    %% Physical location of the cells, assuming located on the same layer, thus
    % spreading the 2-D space
    % locations of excitatory cells
    rng(2024);
    xE = -Ntwk.Scale/2 + gpuArray.rand(Ntwk.Exct.N,1)*Ntwk.Scale;
    [xE, ~] = sort(xE);
    yE = gpuArray.rand(Ntwk.Exct.N,1)*100 - 50; % zeros(Ntwk.Exct.N,1); % assume all neurons approximately on the same layer
    Ntwk.Exct.Location = [xE, yE];
    clear xE yE;
    % locations of Inhibitory cells
    rng(2023);
    xI = -Ntwk.Scale/2 + gpuArray.rand(Ntwk.Inhbt.N,1)*Ntwk.Scale;
    [xI, ~] = sort(xI);
    yI = gpuArray.rand(Ntwk.Inhbt.N,1)*80  - 40;
    Ntwk.Inhbt.Location = [xI, yI];
    clear xI yI;
    if show
        h = figure;
        filename = 'Network structure';
        hold on;
        plot(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 'kv', 'MarkerSize', 4, 'MarkerFaceColor','auto');
        plot(Ntwk.Inhbt.Location(:,1), Ntwk.Inhbt.Location(:,2), '.', 'Color', [.5, .5, .5], 'MarkerSize', 4);
        xlabel('\mum');
        ylabel('\mum');
        xlim([-Ntwk.Scale/2, Ntwk.Scale/2]);
        ylim([-50, 50]);
    end
    % Establish possible synaptic connections, independent from synaptic weights
    % possible synaptic connection from E -> I
    [XE, XI] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1)); % rows represent Inhbt and columns represent Exct
    [YE, YI] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Inhbt.Location(:,2));
    DstcEI = sqrt((XE - XI).^2 + (YE - YI).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
    p_EI = exp(-.5*(DstcEI/Ntwk.Exct.AxonRange).^2); %1 - (1 - exp(-DstcEI/Ntwk.Exct.AxonRange)).^4; % probability of physical connection from E to I based on distance
    Ntwk.Cnnct_EI = p_EI >= rand(size(p_EI)); % connections from E to I, 0 or 1
    clear XE XI YE YI DstcEI p_EI;
    % possible synaptic connection from I -> E
    [XI, XE] = meshgrid(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent Exct and columns represent Inhbt
    [YI, YE] = meshgrid(Ntwk.Inhbt.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcIE = sqrt((XI - XE).^2 + (YI - YE).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
    p_IE = exp(-.5*(DstcIE/Ntwk.Inhbt.AxonRange).^2); % probability of connection from I to E based on distance
    Ntwk.Cnnct_IE = p_IE >= rand(size(p_IE)); % connections from I to E, 0 or 1
    clear XE XI YE YI DstcIE p_IE;
    % possible synaptic connection from E -> E
    [XE1, XE2] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1));
    [YE1, YE2] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcEE = sqrt((XE1 - XE2).^2 + (YE1 - YE2).^2); %  Euclidean distance between each pair of Exct neurons
    p_EE = exp(-.5*(DstcEE/Ntwk.Exct.AxonRange).^2); % probability of connection from E to E based on distance
    Ntwk.Cnnct_EE = p_EE >= rand(size(p_EE)); % connections from E to E, 0 or 1
    clear XE1 XE2 YE1 YE2 DstcEE p_EE;
    % possible synaptic connection from I -> I was not assumed since we do
    % not assume disinhibition in the current interneuronal type (PV like neurons)
    if show
        exampleE = Ntwk.Exct.N/2; %280;
        EIs = find(Ntwk.Cnnct_EI(:, exampleE));
        for i = EIs'
            lgd1 = plot([Ntwk.Exct.Location(exampleE,1), Ntwk.Inhbt.Location(i,1)],...
                [Ntwk.Exct.Location(exampleE,2), Ntwk.Inhbt.Location(i,2)],...
                '-', 'Color', OKeeffe(3,:), 'LineWidth',1);
        end
        %plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        exampleI = Ntwk.Inhbt.N/2;% 48;
        IEs = find(Ntwk.Cnnct_IE(:, exampleI));
        for i = IEs'
            lgd2 = plot([Ntwk.Exct.Location(i,1), Ntwk.Inhbt.Location(exampleI,1)],...
                [Ntwk.Exct.Location(i,2), Ntwk.Inhbt.Location(exampleI,2)],...
                '-', 'Color', OKeeffe(4,:), 'LineWidth',1);
        end
        % plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        legend([lgd1, lgd2], {'E to I','I to E'});
        mysavefig(h, filename, gnrloutdir, 14, [5, 2.5], 2);
        clear EIs IEs i;
    end
    %% Synaptic weights for those possible connections defined above
    % initial weights
    rng(2024);
    Ntwk.wEI_initial = .1*gpuArray.rand(size(Ntwk.Cnnct_EI)).*Ntwk.Cnnct_EI; % synaptic weight from E to I, weak initial connections from E to I
    Ntwk.wIE_initial = .1*gpuArray.rand(size(Ntwk.Cnnct_IE)).*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST
    Ntwk.wEE_initial = .1*gpuArray.rand(size(Ntwk.Cnnct_EE)).*Ntwk.Cnnct_EE; % synaptic weight from E to E, weak initial connections of self-excitation
    
    if show
        h = figure;
        filename = 'wEE_initial';
        imagesc(Ntwk.wEE_initial);
        caxis([0, 1]);
        colormap(bluewhitered(256));
        c = colorbar;
        c.Label.String = 'Initial synaptic weights';
        c.Location = 'northoutside';
        xlabel('Exct channels');
        ylabel('Exct channels');
        mysavefig(h, filename, gnrloutdir, 12, [3.2, 3.6]); % [4, 3]

        h = figure;
        filename = 'wEI_initial';
        imagesc(Ntwk.wEI_initial);
        caxis([0, 1]);
        colormap(bluewhitered(256));
        xlabel('Exct channels');
        ylabel('Inhbt channels');
        mysavefig(h, filename, gnrloutdir, 12, [3.2, 1.46]);

        h = figure;
        filename = 'wIE_initial';
        imagesc(Ntwk.wIE_initial);
        caxis([0, 1]);
        colormap(bluewhitered(256));
        ylabel('Exct channels');
        xlabel('Inhbt channels');
        mysavefig(h, filename, gnrloutdir, 12, [1.6, 3]);
    end
    save(Ntwkfile, 'Ntwk', 'exampleE', 'exampleI');
else
    load(Ntwkfile);
end
