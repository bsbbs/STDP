function InputTuning1D(Ntwk, OKeeffe, plotdir)
%% Input tuning - theoritical
h = figure;
filename = 'InputTuning_Theoritical';
xvec = -Ntwk.XScale:Ntwk.XScale;
subplot(2,1,1); hold on;
lgd = [];
legendLabels = {};
for i = 1:Ntwk.Input.Source
    DstcInput = abs(Ntwk.Input.XLoc(i) - xvec); % Euclidean distance between each pair of neurons
    p_Input = Ntwk.CnnctProb.Input*exp(-.5*(DstcInput/Ntwk.AxonRange.Input).^2); % probability of physical connection based on distance
    y = p_Input;
    y2 = zeros(size(y));
    lgd(i) = plot(xvec, y, 'LineWidth', 2, 'Color', OKeeffe(i,:));
    fill([xvec fliplr(xvec)], [y fliplr(y2)], OKeeffe(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    legendLabels{i} = sprintf('Input %d', i);
end
legend(lgd, legendLabels, 'Location', 'best');
xlabel('Location x (\mum)');
ylabel('Coupling %');
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

subplot(2,1,2); hold on;
undersample = 8;
plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'k^', 'MarkerSize', 3);
plot(Ntwk.Inhbt.Location(1:undersample:end,1), Ntwk.Inhbt.Location(1:undersample:end,2), 'r.', 'MarkerSize', 5);
xlabel('Location x (\mum)');
ylabel('y (\mum)');
xlim([-Ntwk.XScale, Ntwk.XScale]);
ylim([-Ntwk.YScale, Ntwk.YScale]);
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);
%% Input tuning - Empirical 
xinterval = Ntwk.XScale/120;
xwindow = Ntwk.XScale/40;
xvec = -Ntwk.XScale:xinterval:Ntwk.XScale;
Tuning.Input = zeros(Ntwk.Input.Source, numel(xvec));
for xi = 1:numel(xvec)
    for si = 1:Ntwk.Input.Source
        Tuning.Input(si, xi) = sum(Ntwk.Cnnct_Input(Ntwk.Exct.Location(:,1) >= xvec(xi) - xwindow &...
            Ntwk.Exct.Location(:,1) < xvec(xi) + xwindow, Ntwk.Input.Origins == si), 'all');
    end
end

h = figure;
filename = 'InputTuning_Empirical';
subplot(2,1,1); hold on;
lgd = [];
legendLabels = {};
for i = 1:Ntwk.Input.Source
    y = Tuning.Input(i,:)/max(Tuning.Input(:));
    y2 = zeros(size(y));
    lgd(i) = plot(xvec, y, 'LineWidth', 2, 'Color', OKeeffe(i,:));
    fill([xvec fliplr(xvec)], [y fliplr(y2)], OKeeffe(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    legendLabels{i} = sprintf('Input %d', i);
end
legend(lgd, legendLabels, 'Location', 'best');
xlabel('Location x (\mum)');
ylabel('Coupling (normalized)');
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

subplot(2,1,2); hold on;
undersample = 8;
plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'k^', 'MarkerSize', 3);
plot(Ntwk.Inhbt.Location(1:undersample:end,1), Ntwk.Inhbt.Location(1:undersample:end,2), 'r.', 'MarkerSize', 5);
xlabel('Location x (\mum)');
ylabel('y (\mum)');
xlim([-Ntwk.XScale, Ntwk.XScale]);
ylim([-Ntwk.YScale, Ntwk.YScale]);
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

%
filename = 'InputTuningwithExmplNeurons';
ECnnct1 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 1),2);
ECnnct2 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 2),2);
E1E2 = sum(Ntwk.Cnnct_EE(:,ECnnct1 > 0), 2); % Available connections from E1 to E2
E2E1 = sum(Ntwk.Cnnct_EE(ECnnct2 > 0, :)', 2); % Available connections from E2 to E1
E1 = find(ECnnct1(E1E2 & E2E1) > 0 & ECnnct2(E1E2 & E2E1) == 0);
tmp = ECnnct1(E1);
E1 = E1(tmp == max(tmp)); % E1: receiving Input 1 exclusively, but reciprocal connections between E1 and E2
ICnnct1 = sum(Ntwk.Cnnct_EI(:, E1), 2);
ICnnct2 = sum(Ntwk.Cnnct_EI(:, ECnnct2 > 0), 2);
Cnnct1I = sum(Ntwk.Cnnct_IE(E1, :)', 2);
Cnnct2I = sum(Ntwk.Cnnct_IE(ECnnct2 > 0, :)', 2);
I1 = find(ICnnct1 > 0 & ICnnct2 > 0 & Cnnct1I > 0 & Cnnct2I > 0);
tmp = (sum(Ntwk.wEI_initial(I1,ECnnct1>0), 2)).^2; % + (sum(Ntwk.wIE_initial(E1, I1)', 2)).^2;
I1 = I1(tmp == max(tmp)); % I1: receiving from and inhibiting E1 directly, but also tuning to and inhibiting input 2
E2 = find(ECnnct2(E1E2 & E2E1) > 0 & ECnnct1(E1E2 & E2E1) == 0);
tmp = ECnnct2(E2);
E2 = E2(tmp == max(tmp)); % E2: receiving Input 2 exclusively, but reciprocal connections between E1 and E2
ICnnct1 = sum(Ntwk.Cnnct_EI(:, ECnnct1 > 0), 2);
ICnnct2 = sum(Ntwk.Cnnct_EI(:, E2), 2);
Cnnct1I = sum(Ntwk.Cnnct_IE(ECnnct1 > 0, :)', 2);
Cnnct2I = sum(Ntwk.Cnnct_IE(E2, :)', 2);
I2 = find(ICnnct1 > 0 & ICnnct2 > 0 & Cnnct1I > 0 & Cnnct2I > 0);
tmp = (sum(Ntwk.wEI_initial(I2, ECnnct2>0), 2) ).^2; % + (sum(Ntwk.wIE_initial(E2, I2)', 2)).^2;
I2 = I2(tmp == max(tmp)); % I2: receiving from and inhibiting E2 directly, but also tuning to and inhibiting input 1

ICnnct1 = sum(Ntwk.Cnnct_EI(:, ECnnct1 > 0), 2);
ICnnct2 = sum(Ntwk.Cnnct_EI(:, ECnnct2 > 0), 2);
Cnnct1I = sum(Ntwk.Cnnct_IE(ECnnct1 > 0, :)', 2);
Cnnct2I = sum(Ntwk.Cnnct_IE(ECnnct2 > 0, :)', 2);
IShare = find(ICnnct1 > 0 & ICnnct2 > 0 & Cnnct1I > 0 & Cnnct2I > 0);
tmp =  (sum(Ntwk.Cnnct_EI(IShare, ECnnct1 > 0), 2)).^2 + (sum(Ntwk.Cnnct_IE(ECnnct2 > 0, IShare)', 2)).^2;
IShare = IShare(tmp == max(tmp)); % IShare receives projections from both clusters 1 and 2, and inhibiting both clusters
EShare1 = find(Ntwk.Cnnct_EI(IShare,:)' & ECnnct1 > 0);
tmp = Ntwk.wEI_initial(IShare, EShare1)'/max(Ntwk.wEI_initial(IShare, EShare1)) ...
    + sum(Ntwk.Cnnct_Input(EShare1, :), 2)/max(sum(Ntwk.Cnnct_Input(EShare1, :), 2));
EShare1 = EShare1(tmp == max(tmp)); % EShare1 project to IShare and receives from input1
EShare2 = find(Ntwk.Cnnct_IE(:,IShare) & ECnnct2 > 0);
tmp = Ntwk.wIE_initial(EShare2, IShare)/max(Ntwk.wIE_initial(EShare2, IShare)) + sum(Ntwk.Cnnct_Input(EShare2, :), 2)/max(sum(Ntwk.Cnnct_Input(EShare2, :), 2));
EShare2 = EShare2(tmp == max(tmp)); % EShare1 receives inhibition from IShare and receives from input2

Ntwk.Smpl.E = [E1(1), E2(1), EShare1(1), EShare2(1)];
Ntwk.Smpl.I = [I1(1), I2(1), IShare(1)];
for ei = 1:2
    plot(Ntwk.Exct.Location(Ntwk.Smpl.E(ei),1), Ntwk.Exct.Location(Ntwk.Smpl.E(ei),2), 'k^', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end
for ii = 1:2
    plot(Ntwk.Inhbt.Location(Ntwk.Smpl.I(ii),1), Ntwk.Inhbt.Location(Ntwk.Smpl.I(ii),2), 'r.', 'MarkerSize', 14);
end
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

%% %% Visualize connection and pools of samples
Sum1 = sum(Ntwk.wInput(:, Ntwk.Input.Origins == 1), 2);
Sum2 = sum(Ntwk.wInput(:, Ntwk.Input.Origins == 2), 2);
Cnnct1 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 1), 2);
Cnnct2 = sum(Ntwk.Cnnct_Input(:, Ntwk.Input.Origins == 2), 2);
h = figure;
filename = 'Inputtuning_individualE';
plot(Sum1, Sum2, 'k.', 'MarkerSize', 4);
xlim([0,30]);
ylim([0,30]);
xlabel('Tuning weight to Input 1');
ylabel('Tuning weight to Input 2');
mysavefig(h, filename, plotdir, 12, [2.5, 2.2161], 1);

Mat = [sum(~Cnnct1 & ~Cnnct2), sum(Cnnct1 & ~Cnnct2);
    sum(Cnnct2 & ~Cnnct1), sum(Cnnct1 & Cnnct2)];
h = figure;
filename = 'InputTuningNneurons';
imagesc(Mat);
set(gca, 'YDir', 'normal');
colormap('gray');
for i = 1:2
    for j = 1:2
    text(i,j, sprintf('%i',Mat(j,i)),'Color','r');
    end
end
ylabel('Tuning to Input 2');
xlabel('Tuning to Input 1');
mysavefig(h, filename, plotdir, 12, [2.5, 2.2161]); 

h = figure; hold on;
filename = 'InputTuning2D';
scatter(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 3, Sum1, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(Ntwk.Exct.Location(:,1), Ntwk.Exct.Location(:,2), 3, -Sum2, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Tuning weights';
c.Location = 'northoutside';
xlabel('Location x (\mum)');
ylabel('y (\mum)');
mysavefig(h, filename, plotdir, 12, [5, 3], 0);

h = figure; hold on;
filename = 'InputTuning2DSlcted';
TuneMask1 = Cnnct1 & ~Cnnct2;
TuneMask2 = ~Cnnct1 & Cnnct2;
scatter(Ntwk.Exct.Location(TuneMask1,1), Ntwk.Exct.Location(TuneMask1,2), 3, Sum1(TuneMask1), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
scatter(Ntwk.Exct.Location(TuneMask2,1), Ntwk.Exct.Location(TuneMask2,2), 3, -Sum2(TuneMask2), 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
colormap(bluewhitered(256));
c = colorbar;
c.Label.String = 'Tuning weights';
c.Location = 'northoutside';
xlabel('Location x (\mum)');
ylabel('y (\mum)');
mysavefig(h, filename, plotdir, 12, [5, 3], 0);
%% E and I neurons demo
h = figure;
filename = 'EandIPools';
subplot(1,2,1);
plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'k^', 'MarkerSize', 3);
xlabel('Location x (\mum)');
ylabel('y (\mum)');
xlim([-Ntwk.XScale, Ntwk.XScale]);
ylim([-Ntwk.YScale, Ntwk.YScale]);
mysavefig(h, filename, plotdir, 12, [5, 2.5], 2);
subplot(1,2,2);
plot(Ntwk.Inhbt.Location(1:undersample:end,1), Ntwk.Inhbt.Location(1:undersample:end,2), 'r.', 'MarkerSize', 5);
xlabel('Location x (\mum)');
ylabel('y (\mum)');
xlim([-Ntwk.XScale, Ntwk.XScale]);
ylim([-Ntwk.YScale, Ntwk.YScale]);
mysavefig(h, filename, plotdir, 12, [5, 2.5], 2);

