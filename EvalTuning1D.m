function EvalTuning1D(Ntwk, WEE, WEI, WIE, OKeeffe, plotdir)
% Evaluating connection weigths over tuning of inputs.
%% general parameters
xinterval = Ntwk.XScale/60;
xvec = -Ntwk.XScale:xinterval:Ntwk.XScale;
%% Connection weights as a function of the location of the neurons 
% E-to-I connection weights
h = figure; 
filename = 'SynapticWeights_over_Locations';
subplot(3,1,1);hold on;
legendLabels = {'Before','After'};
lgd = [];
for i = 1:2
    if i == 1
        wEI = Ntwk.wEI_initial;
    else
        wEI = WEI;
    end
    Cnnct.EtoI = zeros(1, numel(xvec)-1);
    for xi = 1:(numel(xvec)-1)
        Cnnct.EtoI(xi) = sum(wEI(Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
            Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval, :), 'all');
    end
    if i == 1
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoI, ':', 'LineWidth', .75,'Color', 'k');
    else
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoI, '-', 'LineWidth', .75,'Color', 'k');
    end
end
legend(flip(lgd), flip(legendLabels), 'Location', 'best');
xlabel('Locations of I neurons (\mum)');
ylabel('Effective weights');
title('E to I');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);

% I-to-E connection weights
subplot(3,1,2);hold on;
for i = 1:2
    if i == 1
        wIE = Ntwk.wIE_initial;
    else
        wIE = WIE;
    end
    Cnnct.ItoE = zeros(1, numel(xvec)-1);
    for xi = 1:(numel(xvec)-1)
        Cnnct.ItoE(xi) = sum(wIE(:, Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
            Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval), 'all');
    end
    if i == 1
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.ItoE, ':', 'LineWidth', .75,'Color', 'k');
    else
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.ItoE, '-', 'LineWidth', .75,'Color', 'k');
    end
end
xlabel('Locations of I neurons (\mum)');
ylabel('Synaptic weights');
title('I to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);

% E-to-E connection weights (supplementary figure)
subplot(3,1,3); hold on;
for i = 1:2
    if i == 1
        wEE = Ntwk.wEE_initial;
    else
        wEE = WEE;
    end
    Cnnct.EtoE = zeros(1, numel(xvec)-1);
    for xi = 1:(numel(xvec)-1)
        Cnnct.EtoE(xi) = sum(wEE(:, Ntwk.Exct.Location(:,1) >= xvec(xi) &...
            Ntwk.Exct.Location(:,1) <= xvec(xi) + xinterval), 'all');
    end
    if i == 1
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoE, ':', 'LineWidth', .75, 'Color', 'k');
    else
        lgd(i) = plot(xvec(1:end-1)+xinterval/2, Cnnct.EtoE, '-', 'LineWidth', .75, 'Color', 'k');
    end
end
xlabel('Locations of E neurons (\mum)');
ylabel('Synaptic weights');
title('E to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);

%% Effective connection weights
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
    Tuning.EtoI = zeros(Ntwk.Input.Source, numel(xvec)-1);
    for si = 1:Ntwk.Input.Source
        for xi = 1:(numel(xvec)-1)
            InputtoI = wEI*Ntwk.Cnnct_Input; % from
            Tuning.EtoI(si, xi) = sum(InputtoI(Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
                Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval, Ntwk.Input.Origins == si), 'all');
        end
        if i == 1
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoI(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoI(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of I neurons (\mum)');
ylabel('Effective weights');
title('E to I');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
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
    for si = 1:Ntwk.Input.Source
        for xi = 1:(numel(xvec)-1)
            ItoInput = Ntwk.Cnnct_Input'*wIE; % to
            Tuning.ItoE(si, xi) = sum(ItoInput(Ntwk.Input.Origins == si, Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
                Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval), 'all');
        end
        if i == 1
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.ItoE(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.ItoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of I neurons (\mum)');
ylabel('Effective weights');
title('I to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
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
    for si = 1:Ntwk.Input.Source
        for xi = 1:(numel(xvec)-1)
            EtoInput = Ntwk.Cnnct_Input'*wEE; % to
            Tuning.EtoE(si, xi) = sum(EtoInput(Ntwk.Input.Origins == si, Ntwk.Exct.Location(:,1) >= xvec(xi) &...
                Ntwk.Exct.Location(:,1) <= xvec(xi) + xinterval), 'all');
        end
        if i == 1
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoE(si,:),':', 'LineWidth', .75,'Color', OKeeffe(si,:).^6);
        else
            lgd(i,si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
        end
    end
end
xlabel('Locations of E neurons (\mum)');
ylabel('Effective weights');
title('E to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 2);
%% Effective connection probability
% E-to-I connections, tuning to input sources
h = figure;
filename = 'Effctv_Connections';
subplot(3,1,1); hold on;
legendLabels = {'Input 1','Input 2'};
wEI = Ntwk.Cnnct_EI;
Tuning.EtoI = zeros(Ntwk.Input.Source, numel(xvec)-1);
lgd = [];
for si = 1:Ntwk.Input.Source
    for xi = 1:(numel(xvec)-1)
        InputtoI = wEI*Ntwk.Cnnct_Input; % from
        Tuning.EtoI(si, xi) = sum(InputtoI(Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
            Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval, Ntwk.Input.Origins == si), 'all');
    end
    lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoI(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
end

xlabel('Locations of I neurons (\mum)');
ylabel('Effective connections');
title('E to I');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
legend(lgd, legendLabels, 'Location', 'best');
mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);
% I to E connections
subplot(3,1,2); hold on;
wIE = Ntwk.Cnnct_IE;
Tuning.ItoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
lgd = [];
for si = 1:Ntwk.Input.Source
    for xi = 1:(numel(xvec)-1)
        ItoInput = Ntwk.Cnnct_Input'*wIE; % to
        Tuning.ItoE(si, xi) = sum(ItoInput(Ntwk.Input.Origins == si, Ntwk.Inhbt.Location(:,1) >= xvec(xi) &...
            Ntwk.Inhbt.Location(:,1) <= xvec(xi) + xinterval), 'all');
    end
    lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.ItoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
end
xlabel('Locations of I neurons (\mum)');
ylabel('Effective connections');
title('I to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);
% E to E connections
subplot(3,1,3); hold on;
wEE = Ntwk.Cnnct_EE;
Tuning.EtoE = zeros(Ntwk.Input.Source, numel(xvec)-1);
lgd = [];
for si = 1:Ntwk.Input.Source
    for xi = 1:(numel(xvec)-1)
        EtoInput = Ntwk.Cnnct_Input'*wEE; % to
        Tuning.EtoE(si, xi) = sum(EtoInput(Ntwk.Input.Origins == si, Ntwk.Exct.Location(:,1) >= xvec(xi) &...
            Ntwk.Exct.Location(:,1) <= xvec(xi) + xinterval), 'all');
    end
    lgd(si) = plot(xvec(1:end-1)+xinterval/2, Tuning.EtoE(si,:), '-', 'LineWidth', .75,'Color', OKeeffe(si,:));
end
xlabel('Locations of E neurons (\mum)');
ylabel('Effective connections');
title('E to E');
autoy = ylim();
ylim([autoy(1)-0.1*(range(autoy)),autoy(2)]);
mysavefig(h, filename, plotdir, 12, [2.5, 5], 1);
%% Summed tuning of the feedback inhibition, E_to_I * I_to_E

% initial feedback gain control
InputtoI = Ntwk.wEI_initial*Ntwk.Cnnct_Input;
ItoInput = Ntwk.Cnnct_Input'*Ntwk.wIE_initial;
RecurrentInhbt = ItoInput*InputtoI;
RIMtrx1 = zeros(Ntwk.Input.Source, Ntwk.Input.Source);
for i = 1:2
    for j = 1:2
        RIMtrx1(i,j) = sum(RecurrentInhbt(Ntwk.Input.Origins == i, Ntwk.Input.Origins == j), 'all');
    end
end
% trained feedback gain control
InputtoI = WEI*Ntwk.Cnnct_Input;
ItoInput = Ntwk.Cnnct_Input'*WIE;
RecurrentInhbt = ItoInput*InputtoI;
RIMtrx2 = zeros(Ntwk.Input.Source, Ntwk.Input.Source);
for i = 1:2 % to
    for j = 1:2 % from
        RIMtrx2(i,j) = sum(RecurrentInhbt(Ntwk.Input.Origins == i, Ntwk.Input.Origins == j), 'all');
    end
end
%%
h = figure;
filename = 'FeedbackInhibition';
subplot(1,2,1);
b1 = bar(RIMtrx1,'FaceColor','flat');
ylim([0, max([RIMtrx2(:); RIMtrx1(:)])*1.2]);
xlabel('Origin');
ylabel('Gain control');
ld = legend({'to 1','to 2'}, "Location","best");
mysavefig(h, filename, plotdir, 12, [2.5*2, 2.2161], 2);
subplot(1,2,2);
b2 = bar(RIMtrx2,'FaceColor','flat');
ylim([0, max([RIMtrx2(:); RIMtrx1(:)])*1.2]);
xlabel('Origin');
ylabel('Gain control');
mysavefig(h, filename, plotdir, 12, [2.5*2, 2.2161], 2);





