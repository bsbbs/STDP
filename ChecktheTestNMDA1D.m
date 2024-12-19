% Seqfile = fullfile(plotdir, 'Seq.mat');
% load(Seqfile);
% Ntwkfile = fullfile(plotdir, 'Ntwk.mat');
% load(Ntwkfile);
% Rsltfile = fullfile(subplotdir,'Rslts.mat');
% load(Rsltfile);
%% visualizing dynamically changing periods
smplonsets = round(evs(round(linspace(1,numel(evs(:,1)),30)),1)*1000/dt);
EIMtrx = nan(numel(smplonsets)+1, Ntwk.Input.Source);
IEMtrx = nan(numel(smplonsets)+1,  Ntwk.Input.Source);
RIMtrx = nan(numel(smplonsets)+1, Ntwk.Input.Source, Ntwk.Input.Source);
for evi = 0:numel(smplonsets)
    if evi == 0
        WEI = Ntwk.wEI_initial;
        WIE = Ntwk.wIE_initial;
        WEE = Ntwk.wEE_initial;
    else
        filename = sprintf('RealtimeMonitor_Event%i', evi);
        load(fullfile(subplotdir, [filename, '.mat']));
    end
    InputtoI = WEI*Ntwk.Cnnct_Input;
    ItoInput = Ntwk.Cnnct_Input'*WIE;
    RecurrentInhbt = ItoInput*InputtoI;

    for i = 1:Ntwk.Input.Source
        EIMtrx(evi+1,i) = mean(sum(InputtoI(:,Ntwk.Input.Origins == i), 1));
        IEMtrx(evi+1,i) = mean(sum(ItoInput(Ntwk.Input.Origins == i,:), 2));
    end
    for i = 1:Ntwk.Input.Source % to
        for j = 1:Ntwk.Input.Source % from
            RIMtrx(evi+1,i,j) = mean(RecurrentInhbt(Ntwk.Input.Origins == i, Ntwk.Input.Origins == j), 'all');
        end
    end
end

h = figure;
filename = 'WeightsSummaryDynamics';
subplot(2,2,1); hold on;
for i = 1:2
    plot(0:18,EIMtrx(:,i), 'Color',OKeeffe(i,:));
end
legend({'E to I from Input 1','E to I from Input 2'}, 'Location','best');
xlabel('Check time point');
ylabel('Summed weights');
mysavefig(h, filename, subplotdir, 12, [5, 5], 1);
subplot(2,2,2); hold on;
for i = 1:2
    plot(0:18,IEMtrx(:,i), 'Color',OKeeffe(i,:));
end
legend({'I to E to Input 1','I to E to Input 2'}, 'Location','best');
xlabel('Check time point');
ylabel('Summed weights');
mysavefig(h, filename, subplotdir, 12, [5, 5], 1);
subplot(2,2,3); hold on;
for i = 1:2
    plot(0:18,RIMtrx(:,i,1), 'Color',OKeeffe(i,:));
end
legend({'Input 1 to Input 1','Input 1 to Input 2'}, 'Location','best');
xlabel('Check time point');
ylabel('Summed weights');
mysavefig(h, filename, subplotdir, 12, [5, 5], 1);
subplot(2,2,4); hold on;
for i = 1:2
    plot(0:18,RIMtrx(:,i,2), 'Color',OKeeffe(i,:));
end
legend({'Input 2 to Input 1','Input 2 to Input 2'}, 'Location','best');
xlabel('Check time point');
ylabel('Summed weights');
mysavefig(h, filename, subplotdir, 12, [5, 5], 1);
%%
for evi = 1:numel(smplonsets)
    filename = sprintf('RealtimeMonitor_Event%i', evi);
    load(fullfile(subplotdir, [filename, '.mat']));
    plotdirintersect = fullfile(subplotdir, sprintf('Intersect%i', evi));
    if ~exist(plotdirintersect,'dir')
        mkdir(plotdirintersect);
    end
    %%
    h = figure;
    filename = 'Example neurons activity';
%     subplot(2,1,1); hold on;
%     lg = [];
%     lg(1) = plot(timevec - min(timevec), Smpl.ExctgAMPA(1:numel(timevec),1), '-');
%     lg(2) = plot(timevec - min(timevec), Smpl.ExctgNMDA(1:numel(timevec),1), '-');
%     lg(3) = plot(timevec - min(timevec), Smpl.ExctgGABA(1:numel(timevec),1), '-');
%     lg(4) = plot(timevec - min(timevec), Smpl.InhbtgAMPA(1:numel(timevec),1), '-');
%     lg(5) = plot(timevec - min(timevec), Smpl.InhbtgNMDA(1:numel(timevec),1), '-');
%     legend(lg, {'gAMPA on E','gNMDA on E','gGABA on E', 'gAMPA on I','gNMDA on I'}, 'Location','best');
%     xlabel('Time (ms)');
%     ylabel('Conductance (nS)');
%     title('Synaptic activity');
%     mysavefig(h, filename, plotdirintersect, 12, [4,4], 1);
    subplot(4,1,1); hold on;
    ExctP = Smpl.ExctV(1:numel(timevec),1);
    ExctP(ExctP == Ntwk.Vfire) = 0;
    plot(timevec - min(timevec), ExctP, 'k-');
    xlim([95, 600]);
    ylim([-80, -0]);
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,2); hold on;
    InhbtP = Smpl.InhbtV(1:numel(timevec),1);
    InhbtP(InhbtP == Ntwk.Vfire) = 0;
    plot(timevec - min(timevec), InhbtP, 'r-');
    %legend(lg, {'Exct', 'Inhbt'}, 'Location','best');
    xlim([95, 600]);
    ylim([-80, -0]);
    xlabel('Time (ms)');
    ylabel('Membrane potential (mV)');
    mysavefig(h, filename, plotdirintersect, 12, [4, 8], 1);
    subplot(4,1,3); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WEI(1:numel(timevec),1,1)), 'k-');
    xlim([95, 600]);
    ylim([.8847, .8857]);
    xlabel('Time (ms)');
    ylabel('wEI');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,4); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WIE(1:numel(timevec),1,1)), 'r-');
    xlim([95, 600]);
    ylim([.976, .982]);
    xlabel('Time (ms)');
    ylabel('wIE');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    %% Synaptic weights on example neurons
    h = figure;
    filename = 'Example synaptic weights';
    subplot(4,1,1); hold on;
    lg = [];
    lg(1) = plot(timevec - min(timevec), Smpl.xEpre(1:numel(timevec),1), 'k-');
    lg(2) = plot(timevec - min(timevec), Smpl.xIpre(1:numel(timevec),1), 'r-');
    xlim([95, 200]);
    xlabel('Time (ms)');
    ylabel('Integration');
    title('STDP Integration');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,2); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WEI(1:numel(timevec),1,1)), 'k-');
    xlim([95, 200]);
    xlabel('Time (ms)');
    ylabel('wEI');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,3); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WIE(1:numel(timevec),1,1)), 'r-');
    xlim([95, 200]);
    xlabel('Time (ms)');
    ylabel('wIE');
    mysavefig(h, filename, plotdirintersect, 12, [4,8], 1);
    subplot(4,1,4); hold on;
    plot(timevec - min(timevec), squeeze(Smpl.WEE(1:numel(timevec),2,1)), 'k--');
    xlim([95, 200]);
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
%% old version
%% visualizing example neurons
% % load(Rsltfile);
% h = figure;
% filename = 'Example neurons activity';
% subplot(2,1,1); hold on;
% lg = [];
% lg(1) = plot(time/1000, Smpl.ExctgE(:,1), 'k-');
% lg(2) = plot(time/1000, Smpl.ExctgI(:,1), 'k--');
% lg(3) = plot(time/1000, Smpl.InhbtgE(:,1), 'r-');
% lg(4) = plot(time/1000, Smpl.InhbtgI(:,1), 'r--');
% legend(lg, {'gE on Exct cell','gI on Exct cell', 'gE on Inhbt cell','gI on Inhbt cell'}, 'Location','best');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('Conductance (nS)');
% title('Synaptic activity');
% mysavefig(h, filename, plotdir, 12, [8,4], 1);
% subplot(2,2,3); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.ExctgE(:,1), 'k-');
% lg(2) = plot(time, Smpl.ExctgI(:,1), 'k--');
% lg(3) = plot(time, Smpl.InhbtgE(:,1), 'r-');
% lg(4) = plot(time, Smpl.InhbtgI(:,1), 'r--');
% xlim([leftt, leftt+700]);
% xlabel('Time (ms)');
% ylabel('Conductance (nS)');
% title('Synaptic activity');
% mysavefig(h, filename, plotdir, 12, [8,4], 1);
% subplot(2,2,4); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.ExctV(:,1), 'k-');
% lg(2) = plot(time, Smpl.InhbtV(:,1), 'r-');
% legend(lg, {'Exct', 'Inhbt'}, 'Location','best');
% xlim([leftt, leftt+700]);
% xlabel('Time (ms)');
% ylabel('Membrane potential (mV)');
% title('Activity of a single neuron');
% mysavefig(h, filename, plotdir, 12, [8, 4], 1);
% %% Synaptic weights on example neurons
% h = figure;
% filename = 'Example synaptic weights';
% subplot(4,2,1); hold on;
% lg = [];
% lg(1) = plot(time/1000, Smpl.xEpre(:,1), 'k-');
% lg(2) = plot(time/1000, Smpl.xIpre(:,1), 'r-');
% %lg(3) = plot(time/1000, Exmpl.xEpost, 'k--');
% %lg(4) = plot(time/1000, Exmpl.xIpost, 'r--');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('Integration');
% title('STDP Integration');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,3); hold on;
% plot(time/1000, squeeze(Smpl.WEI(:,1,1)), 'k-');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wEI');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,5); hold on;
% plot(time/1000, squeeze(Smpl.WIE(:,1,1)), 'r-');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wIE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,7); hold on;
% plot(time/1000, squeeze(Smpl.WEE(:,2,1)), 'k--');
% xlim([0 duration/1000]);
% xlabel('Time (s)');
% ylabel('wEE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,2); hold on;
% lg = [];
% lg(1) = plot(time, Smpl.xEpre(:,1), 'k-');
% lg(2) = plot(time, Smpl.xIpre(:,1), 'r-');
% % lg(3) = plot(time, Exmpl.xEpost(1,:), 'k--');
% % lg(4) = plot(time, Exmpl.xIpost(1,:), 'r--');
% legend(lg, {'xE', 'xI'}, 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('Integration');
% title('STDP Integration');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,4); hold on;
% plot(time, squeeze(Smpl.WEI(:,1,1)), 'k-');
% legend('WEI', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wEI');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,6); hold on;
% plot(time, squeeze(Smpl.WIE(:,1,1)), 'r-');
% legend('WIE', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wIE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
% subplot(4,2,8); hold on;
% plot(time, squeeze(Smpl.WEE(:,2,1)), 'k--');
% legend('WEE', 'Location','best');
% xlim([leftt leftt+700]);
% xlabel('Time (ms)');
% ylabel('wEE');
% mysavefig(h, filename, plotdir, 12, [8,8], 1);
    %%
    EvalTuning1D(Ntwk,WEE,WEI,WIE,OKeeffe,plotdirintersect);
    close all;
end

