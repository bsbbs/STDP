function EvalTuning1D(Ntwk, WEE, WEI, WIE, OKeeffe, plotdir)
% Evaluating connection weigths over tuning of inputs.
%% general parameters
xinterval = Ntwk.XScale/60;
xvec = -Ntwk.XScale:xinterval:Ntwk.XScale;
%% target the tuning neurons
Cnnct1 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 1), 2);
Cnnct2 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 2), 2);
Cnnct = sum(Ntwk.Cnnct_Input, 2);
InhbtCnnct = Ntwk.Cnnct_EI*Cnnct;
TuneMasks(:,1) = Cnnct1 & ~Cnnct2;
TuneMasks(:,2) = ~Cnnct1 & Cnnct2;
TuneMask = Cnnct > 3;
TuneMaskInhbt = InhbtCnnct > 200;
%% Connection weights as a function of the location of the neurons 
% % E-to-I connection weights
% h = figure; 
% filename = 'SynapticWeights_over_Locations';
% subplot(3,1,1);hold on;
% legendLabels = {'Before','After'};
% lgd = [];
% for i = 1:2
%     if i == 1
%         wEI = Ntwk.wEI_initial;
%     else
%         wEI = WEI;
%     end
%     Cnnct.EtoI = zeros(1, numel(xvec)-1);
%     for xi = 1:(numel(xvec)-1)
%         Cnnct.EtoI(xi) = sum(wEI(Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
%             Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval, :), 'all');
%     end
%     if i == 1
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoI, ':', 'LineWidth', .75,'Color', 'k');
%     else
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoI, '-', 'LineWidth', .75,'Color', 'k');
%     end
% end
% legend(flip(lgd), flip(legendLabels), 'Location', 'best');
% xlabel('Locations of I neurons (\mum)');
% ylabel('Effective weights');
% title('E to I');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);
% 
% % I-to-E connection weights
% subplot(3,1,2);hold on;
% for i = 1:2
%     if i == 1
%         wIE = Ntwk.wIE_initial;
%     else
%         wIE = WIE;
%     end
%     Cnnct.ItoE = zeros(1, numel(xvec)-1);
%     for xi = 1:(numel(xvec)-1)
%         Cnnct.ItoE(xi) = sum(wIE(:, Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
%             Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval), 'all');
%     end
%     if i == 1
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.ItoE, ':', 'LineWidth', .75,'Color', 'k');
%     else
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.ItoE, '-', 'LineWidth', .75,'Color', 'k');
%     end
% end
% xlabel('Locations of I neurons (\mum)');
% ylabel('Synaptic weights');
% title('I to E');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);
% 
% % E-to-E connection weights (supplementary figure)
% subplot(3,1,3); hold on;
% for i = 1:2
%     if i == 1
%         wEE = Ntwk.wEE_initial;
%     else
%         wEE = WEE;
%     end
%     Cnnct.EtoE = zeros(1, numel(xvec)-1);
%     for xi = 1:(numel(xvec)-1)
%         Cnnct.EtoE(xi) = sum(wEE(:, Ntwk.Exct.Location(:,1) >= xvec(xi) &...
%             Ntwk.Exct.Location(:,1) <= xvec(xi) + xinterval), 'all');
%     end
%     if i == 1
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoE, ':', 'LineWidth', .75, 'Color', 'k');
%     else
%         lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoE, '-', 'LineWidth', .75, 'Color', 'k');
%     end
% end
% xlabel('Locations of E neurons (\mum)');
% ylabel('Synaptic weights');
% title('E to E');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);

%% Effective connection weights
xinterval = Ntwk.XScale/120;
xwindow = Ntwk.XScale/7;
xvec = -Ntwk.XScale:xinterval:Ntwk.XScale;
% E-to-I connection weights, tuning to input sources
h = figure;
filename = 'Effctv_weights';
subplot(3,1,1); hold on;
legendLabels = {'After - Input 1','Before - Input 1', 'After - Input 2','Before - Input 2'};
for i = 1:2
    if i == 1
        wEI = Ntwk.wEI_initial;
    else
        wEI = WEI;
    end
    Tuning.EtoI = zeros(Ntwk.Input.Source, numel(xvec));
    InputtoI = wEI*TuneMasks./sum(TuneMasks, 1); %Ntwk.Cnnct_Input; % from
    for si = 1:Ntwk.Input.Source
        for xi = 1:numel(xvec)
            Tuning.EtoI(si, xi) = mean(InputtoI(Ntwk.Inhbt.Location(:,1) >= xvec(xi)- xwindow &...
                Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xwindow, si));
        end
        if i == 1
            lgd(i,si) = plot(xvec, Tuning.EtoI(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec, Tuning.EtoI(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of I neurons (\mum)');
ylabel('Effective weights');
title('E to I');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% ylim([0, .04]);
tmp = flip(lgd, 1);
legend(tmp(:), legendLabels, 'Location', 'best', 'NumColumns',2);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);
% I to E connections
subplot(3,1,2); hold on;
for i = 1:2
    if i == 1
        wIE = Ntwk.wIE_initial;
    else
        wIE = WIE;
    end
    Tuning.ItoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
    ItoInput = TuneMasks'*wIE./sum(TuneMasks, 1)'; % to
    for si = 1:Ntwk.Input.Source
        for xi = 1:numel(xvec)
            Tuning.ItoE(si, xi) = mean(ItoInput(si, Ntwk.Inhbt.Location(:,1) >= xvec(xi) - xwindow &...
                Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xwindow));
        end
        if i == 1
            lgd(i,si) = plot(xvec, Tuning.ItoE(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec, Tuning.ItoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of I neurons (\mum)');
ylabel('Effective weights');
title('I to E');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% ylim([0, .04]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);
% E to E connections
subplot(3,1,3); hold on;
for i = 1:2
    if i == 1
        wEE = Ntwk.wEE_initial;
    else
        wEE = WEE;
    end
    Tuning.EtoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
    EtoInput = TuneMasks'*wEE./sum(TuneMasks, 1)'; % to
    for si = 1:Ntwk.Input.Source
        for xi = 1:numel(xvec)
            Tuning.EtoE(si, xi) = mean(EtoInput(si, Ntwk.Exct.Location(:,1) >= xvec(xi) - xwindow &...
                Ntwk.Exct.Location(:,1) <= xvec(xi) + xwindow));
        end
        if i == 1
            lgd(i,si) = plot(xvec, Tuning.EtoE(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec, Tuning.EtoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of E neurons (\mum)');
ylabel('Effective weights');
title('E to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);

%% Summed tuning of the feedback inhibition, E_to_I * I_to_E

% initial feedback gain control
InputtoI = Ntwk.wEI_initial*TuneMasks./sum(TuneMasks, 1); % from
ItoInput = TuneMasks'*Ntwk.wIE_initial./sum(TuneMasks, 1)'; % to
RecurrentInhbt1 = (ItoInput*InputtoI);
% trained feedback gain control
InputtoI = WEI*TuneMasks./sum(TuneMasks, 1); % from
ItoInput = TuneMasks'*WIE./sum(TuneMasks, 1)'; % to
RecurrentInhbt2 = ItoInput*InputtoI;
yl = max([RecurrentInhbt2(:); RecurrentInhbt1(:)])*1.2;

h = figure;
filename = 'FeedbackInhibition';
subplot(1,2,2);
b2 = bar(RecurrentInhbt2,'FaceColor','flat');
aftery = ylim;
xlabel('Input origin');
ylabel('Excitation x Gain control');
mysavefig(h, filename, plotdir, 12, [2.5*2, 2.2161], 2);
subplot(1,2,1);
b1 = bar(RecurrentInhbt1,'FaceColor','flat');
ylim(aftery);
xlabel('Input origin');
ylabel('Excitation x Gain control');
ld = legend({'to 1','to 2'}, "Location","best");
mysavefig(h, filename, plotdir, 12, [2.5*2, 2.2161], 2);


%% Effective connection probability
% % E-to-I connections, tuning to input sources
% h = figure;
% filename = 'Effctv_Connections';
% subplot(3,1,1); hold on;
% legendLabels = {'Input 1','Input 2'};
% wEI = Ntwk.Cnnct_EI;
% Tuning.EtoI = zeros(Ntwk.Input.Source, numel(xvec)-1);
% lgd = [];
% for si = 1:Ntwk.Input.Source
%     for xi = 1:(numel(xvec)-1)
%         InputtoI = wEI*Ntwk.Cnnct_Input; % from
%         Tuning.EtoI(si, xi) = sum(InputtoI(Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
%             Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval, Ntwk.Input.Origins == si), 'all');
%     end
%     lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoI(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
% end
% 
% xlabel('Locations of I neurons (\mum)');
% ylabel('Effective connections');
% title('E to I');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% legend(lgd, legendLabels, 'Location', 'best');
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);
% % I to E connections
% subplot(3,1,2); hold on;
% wIE = Ntwk.Cnnct_IE;
% Tuning.ItoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
% lgd = [];
% for si = 1:Ntwk.Input.Source
%     for xi = 1:(numel(xvec)-1)
%         ItoInput = Ntwk.Cnnct_Input'*wIE; % to
%         Tuning.ItoE(si, xi) = sum(ItoInput(Ntwk.Input.Origins == si, Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
%             Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval), 'all');
%     end
%     lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.ItoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
% end
% xlabel('Locations of I neurons (\mum)');
% ylabel('Effective connections');
% title('I to E');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);
% % E to E connections
% subplot(3,1,3); hold on;
% wEE = Ntwk.Cnnct_EE;
% Tuning.EtoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
% lgd = [];
% for si = 1:Ntwk.Input.Source
%     for xi = 1:(numel(xvec)-1)
%         EtoInput = Ntwk.Cnnct_Input'*wEE; % to
%         Tuning.EtoE(si, xi) = sum(EtoInput(Ntwk.Input.Origins == si, Ntwk.Exct.Location(:,1) >= xvec(xi) &...
%             Ntwk.Exct.Location(:,1) <= xvec(xi) + xinterval), 'all');
%     end
%     lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
% end
% xlabel('Locations of E neurons (\mum)');
% ylabel('Effective connections');
% title('E to E');
% autoy = ylim();
% ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
% mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);


