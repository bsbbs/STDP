Ntwkfile = fullfile(gnrloutdir, 'NtwkPrdcGPU.mat');
if ~exist(Ntwkfile, 'file')
    %% The structure of the network
    Ntwk.Scale = 500; % um, scale of the micro structure to test, assumed as one dimentional 
    Ntwk.Exct.Location = []; % the physical location of excitatory cells
    Ntwk.Exct.N = 8000; % number of the excitarory cells
    Ntwk.Exct.AxonRange = 130; 
    Ntwk.Inhbt.Location = [];
    Ntwk.Inhbt.N = 2000;
    Ntwk.Inhbt.AxonRange = 97; % um, the standard deviation of axon physical connection range for interneurons
    Ntwk.AxonRange.EE = 130; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and E
    Ntwk.AxonRange.EI = 100; % um, standard deviation of Gaussian decay of connection probability over somatic distance of E and I 
    Ntwk.AxonRange.IE = 97; % um, standard deviation of Gaussian decay of connection probability over somatic distance of I and E 
    Ntwk.CnnctProb.EE = .1; % the maximum connection probabability from E to E
    Ntwk.CnnctProb.EI = .2; % the maximum connection probabability from E to I
    Ntwk.CnnctProb.IE = .3; % the maximum connection probabability from I to E
    Ntwk.Smpl.E = 3999;
    Ntwk.Smpl.I = 1000;
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
    % rising time of those synapses were assumed as instant
    Ntwk.Synapse.gbarE = .140; % nS, the maximum excitatory synaptic conductance, like AMPA (Vogels et al., 2011)
    Ntwk.Synapse.gbarI = .350; % nS, the maximum inhibitory synaptic conductance, like GABA (Vogels et al., 2011)
    Ntwk.Synapse.gbarE = 4.7; % original value .14 nS
    Ntwk.Synapse.gbarI = 31.5; % original value .35 nS
    %% Spike-timing dependent plasticity
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
    Ntwk.InhbtSTDP.intercept_pre = -2*Ntwk.InhbtSTDP.rho/1000*Ntwk.InhbtSTDP.tau_prepost; % depression if no post spikes
    Ntwk.InhbtSTDP.intercept_post = 0; % depression if no pre spikes
    %% Other parameters
    Ntwk.Noise.tauN = 2; % ms, time constant for the OU process of noise
    Ntwk.Noise.sgm = 30; % amplitude of noise

    %% Physical location of the cells, assuming located on the same layer, thus
    % spreading the 2-D space
    % locations of excitatory cells
    randomseed = 2026;
    gpurng(randomseed);
    xE = -Ntwk.Scale + gpuArray.rand(Ntwk.Exct.N,1)*Ntwk.Scale*2;
    [xE, ~] = sort(xE);
    gpurng(randomseed);
    yE = -Ntwk.Scale/2 + gpuArray.rand(Ntwk.Exct.N,1)*Ntwk.Scale; % zeros(Ntwk.Exct.N,1); % assume all neurons approximately on the same layer
    Ntwk.Exct.Location = [xE, yE];
    clear xE yE;
    % locations of Inhibitory cells
    gpurng(randomseed);
    xI = -Ntwk.Scale + gpuArray.rand(Ntwk.Inhbt.N,1)*Ntwk.Scale*2;
    [xI, ~] = sort(xI);
    gpurng(randomseed);
    yI = -Ntwk.Scale/2 + gpuArray.rand(Ntwk.Inhbt.N,1)*Ntwk.Scale;
    Ntwk.Inhbt.Location = [xI, yI];
    clear xI yI;
    
    % Establish possible synaptic connections, independent from synaptic weights
    % possible synaptic connection from E -> I
    [XE, XI] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Inhbt.Location(:,1)); % rows represent Inhbt and columns represent Exct
    [YE, YI] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Inhbt.Location(:,2));
    DstcEI = sqrt(min(Ntwk.Scale*2-abs(XE - XI), abs(XE - XI)).^2 + min(Ntwk.Scale - abs(YE - YI), abs(YE - YI)).^2); % Euclidean distance between each pair of Exct and Inhbt neurons
    p_EI = Ntwk.CnnctProb.EI*exp(-.5*(DstcEI/Ntwk.AxonRange.EI).^2); %1 - (1 - exp(-DstcEI/Ntwk.Exct.AxonRange)).^4; % probability of physical connection from E to I based on distance
    gpurng(randomseed);
    Ntwk.Cnnct_EI = p_EI >= gpuArray.rand(size(p_EI)); % connections from E to I, 0 or 1
    clear XE XI YE YI DstcEI p_EI;
    % possible synaptic connection from I -> E
    [XI, XE] = meshgrid(Ntwk.Inhbt.Location(:,1), Ntwk.Exct.Location(:,1)); % rows represent Exct and columns represent Inhbt
    [YI, YE] = meshgrid(Ntwk.Inhbt.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcIE = sqrt(min(Ntwk.Scale*2-abs(XI - XE), abs(XI - XE)).^2 + min(Ntwk.Scale - abs(YI - YE), abs(YI - YE)).^2); % Euclidean distance between each pair of Inhbt and Exct neurons
    p_IE = Ntwk.CnnctProb.IE*exp(-.5*(DstcIE/Ntwk.AxonRange.IE).^2); % probability of connection from I to E based on distance
    gpurng(randomseed);
    Ntwk.Cnnct_IE = p_IE >= gpuArray.rand(size(p_IE)); % connections from I to E, 0 or 1
    clear XE XI YE YI DstcIE p_IE;
    % possible synaptic connection from E -> E
    [XE1, XE2] = meshgrid(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,1));
    [YE1, YE2] = meshgrid(Ntwk.Exct.Location(:,2), Ntwk.Exct.Location(:,2));
    DstcEE = sqrt(min(Ntwk.Scale*2-abs(XE1 - XE2), abs(XE1 - XE2)).^2 + min(Ntwk.Scale - abs(YE1 - YE2), abs(YE1 - YE2)).^2); %  Euclidean distance between each pair of Exct neurons
    p_EE = Ntwk.CnnctProb.EE*exp(-.5*(DstcEE/Ntwk.AxonRange.EE).^2); % probability of connection from E to E based on distance
    gpurng(randomseed);
    Ntwk.Cnnct_EE = p_EE >= gpuArray.rand(size(p_EE)); % connections from E to E, 0 or 1
    clear XE1 XE2 YE1 YE2 DstcEE p_EE;
    % possible synaptic connection from I -> I was not assumed since we do
    % not assume disinhibition in the current interneuronal type (PV like neurons)
    %% Visualization
    downsmpl = 48;
    if show
        h = figure;
        filename = 'Network structure';
        hold on;
        plot(Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,1), Ntwk.Exct.Location(1:downsmpl:Ntwk.Exct.N,2), 'k^', 'MarkerSize', 4, 'MarkerFaceColor','auto');
        plot(Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,1),...
            Ntwk.Inhbt.Location(1:round(downsmpl/Ntwk.Exct.N*Ntwk.Inhbt.N):Ntwk.Inhbt.N,2), '.', 'Color', 'r', 'MarkerSize', 4);
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        xlim([-Ntwk.Scale, Ntwk.Scale]);
        ylim([-Ntwk.Scale/2, Ntwk.Scale/2]);
   
        EIs = find(Ntwk.Cnnct_EI(:, Ntwk.Smpl.E));
        for i = EIs'
            lgd1 = plot([Ntwk.Exct.Location(Ntwk.Smpl.E,1), Ntwk.Inhbt.Location(i,1)],...
                [Ntwk.Exct.Location(Ntwk.Smpl.E,2), Ntwk.Inhbt.Location(i,2)],...
                '-', 'Color', 'k', 'LineWidth',1);
        end
        %plot(Ntwk.Inhbt.Location(EIs,1), Ntwk.Inhbt.Location(EIs,2), '.', 'Color', OKeeffe(8,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        
        IEs = find(Ntwk.Cnnct_IE(:, Ntwk.Smpl.I));
        for i = IEs'
            lgd2 = plot([Ntwk.Exct.Location(i,1), Ntwk.Inhbt.Location(Ntwk.Smpl.I,1)],...
                [Ntwk.Exct.Location(i,2), Ntwk.Inhbt.Location(Ntwk.Smpl.I,2)],...
                '-', 'Color', 'r', 'LineWidth',1);
        end
        % plot(Ntwk.Inhbt.Location(IEs,1), Ntwk.Inhbt.Location(IEs,2), '.', 'Color', OKeeffe(7,:), 'MarkerSize', Ntwk.Inhbt.Properties.size);
        legend([lgd1, lgd2], {'E to I','I to E'});
        mysavefig(h, filename, gnrloutdir, 14, [5, 2.5], 2);
        clear EIs IEs i;
    end
    %% Synaptic weights for those possible connections defined above
    % initial weights
    gpurng(randomseed+1);
    Ntwk.wEI_initial = .05*gpuArray.rand(size(Ntwk.Cnnct_EI)).*Ntwk.Cnnct_EI; % synaptic weight from E to I, weak initial connections from E to I
    Ntwk.wIE_initial = .05*gpuArray.rand(size(Ntwk.Cnnct_IE)).*Ntwk.Cnnct_IE; % synaptic weight from I to E, weak initial connections from the nearby SST
    Ntwk.wEE_initial = .05*gpuArray.rand(size(Ntwk.Cnnct_EE)).*Ntwk.Cnnct_EE; % synaptic weight from E to E, weak initial connections of self-excitation
    %% visualization
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
        mysavefig(h, filename, gnrloutdir, 12, [2.5, 2.81]); % [4, 3] % [3.2, 3.6]

        h = figure;
        filename = 'wEI_initial';
        imagesc(Ntwk.wEI_initial);
        caxis([0, 1]);
        colormap(bluewhitered(256));
        xlabel('Exct channels');
        ylabel('Inhbt channels');
        mysavefig(h, filename, gnrloutdir, 12, [2.5, 2.2161]);

        h = figure;
        filename = 'wIE_initial';
        imagesc(Ntwk.wIE_initial);
        caxis([0, 1]);
        colormap(bluewhitered(256));
        ylabel('Exct channels');
        xlabel('Inhbt channels');
        mysavefig(h, filename, gnrloutdir, 12, [2.5, 2.2161]); % 2.5*149.1848/168.2961
    end
    save(Ntwkfile, 'Ntwk');
else
    load(Ntwkfile);
end
