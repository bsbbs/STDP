% Check the test
DefineIO;
Ntrial = 1600;
ProjectName = sprintf('SyncNMDA_%i', Ntrial);
plotdir = fullfile(Projdir, ProjectName);
if ~exist(plotdir,'dir')
    mkdir(plotdir);
end
Seqfile = fullfile(plotdir, 'Seq.mat');
load(Seqfile);
Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
load(Ntwkfile);
%% visualizing dynamically changing periods
smplonsets = round(evs(round(linspace(1,numel(evs(:,1)),18)),1)*1000/dt);
for evi = 1:numel(smplonsets)
    filename = sprintf('RealtimeMonitor_Event%i', evi);
    load(fullfile(plotdir, [filename, '.mat']));
    plotdirintersect = fullfile(plotdir, sprintf('Intersect%i', evi));
    if ~exist(plotdirintersect,'dir')
        mkdir(plotdirintersect);
    end
    %%
    h = figure;
    filename = 'Example neurons activity';
    subplot(2,1,1); hold on;
    lg = [];
    lg(1) = plot(timevec - min(timevec), Smpl.ExctgAMPA(:,1), '-');
    lg(2) = plot(timevec - min(timevec), Smpl.ExctgNMDA(:,1), '-');
    lg(3) = plot(timevec - min(timevec), Smpl.ExctgGABA(:,1), '-');
    lg(4) = plot(timevec - min(timevec), Smpl.InhbtgAMPA(:,1), '-');
    lg(5) = plot(timevec - min(timevec), Smpl.InhbtgNMDA(:,1), '-');
    legend(lg, {'gAMPA on E','gNMDA on E','gGABA on E', 'gAMPA on I','gNMDA on I'}, 'Location','best');
    xlabel('Time (ms)');
    ylabel('Conductance (nS)');
    title('Synaptic activity');
    mysavefig(h, filename, plotdirintersect, 12, [4,4], 1);
    subplot(2,1,2); hold on;
    lg = [];
    lg(2) = plot(timevec - min(timevec), Smpl.InhbtV(:,1), 'r-');
    lg(1) = plot(timevec - min(timevec), Smpl.ExctV(:,1), 'k-');
    legend(lg, {'Exct', 'Inhbt'}, 'Location','best');
    ylim([-70, -45]);
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    title('Activity of a single neuron');
    mysavefig(h, filename, plotdirintersect, 12, [4, 4], 1);
    %% Synaptic weights on example neurons
    h = figure;
    filename = 'Example synaptic weights';
    subplot(4,1,1); hold on;
    lg = [];
    lg(1) = plot(timevec - min(timevec), Smpl.xEpre(:,1), 'k-');
    lg(2) = plot(timevec - min(timevec), Smpl.xIpre(:,1), 'r-');
    xlabel('Time (ms)');
    ylabel('Integration');
    title('STDP Integration');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,2); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WEI(:,1,1)), 'k-');
    xlabel('Time (ms)');
    ylabel('wEI');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,3); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WIE(:,1,1)), 'r-');
    xlabel('Time (ms)');
    ylabel('wIE');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,4); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WEE(:,2,1)), 'k--');
    xlabel('Time (ms)');
    ylabel('wEE');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
  
    %% visualize weight change
    h = figure;
    filename = sprintf('WEI_change_Intersect%i', evi);
    imagesc(WEI - Ntwk.wEI_initial)
    colormap(bluewhitered);
    c = colorbar;
    c.Label.String = 'Weight change';
    c.Location = 'northoutside';
    xlabel("Exct neurons");
    ylabel("Inhbt neurons");
    mysavefig(h, filename, plotdirintersect, 12, [2.5, 2.81], 1);

    h = figure;
    filename = sprintf('WIE_change_Intersect%i', evi);
    imagesc(WIE - Ntwk.wIE_initial)
    colormap(bluewhitered);
    c = colorbar;
    c.Label.String = 'Weight change';
    c.Location = 'northoutside';
    xlabel("Inhbt neurons");
    ylabel("Exct neurons");
    mysavefig(h, filename, plotdirintersect, 12, [2.5, 2.81], 1);

    h = figure;
    filename = sprintf('WEE_change_Intersect%i', evi);
    imagesc(WEE - Ntwk.wEE_initial)
    colormap(bluewhitered);
    c = colorbar;
    c.Label.String = 'Weight change';
    c.Location = 'northoutside';
    xlabel("Exct neurons");
    ylabel("Exct neurons");
    mysavefig(h, filename, plotdirintersect, 12, [2.5, 2.81], 1);
    %%
    EvalTuning(Ntwk,WEE,WEI,WIE,OKeeffe,plotdirintersect);
end