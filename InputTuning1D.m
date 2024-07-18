function Ntwk = InputTuning1D(Ntwk, OKeeffe, plotdir)
%% Input tuning
xinterval = Ntwk.XScale/60;
xvec = -Ntwk.Scale:xinterval:Ntwk.XScale;
Tuning.Input = zeros(Ntwk.Input.Source, numel(xvec));
for xi = 1:numel(xvec)
    for si = 1:Ntwk.Input.Source
        Tuning.Input(si, xi) = sum(Ntwk.Cnnct_Input(Ntwk.Exct.Location(:,1) >= xvec(xi) &...
            Ntwk.Exct.Location(:,1) < xvec(xi) + xinterval, Ntwk.Input.Origins == si), 'all');
    end
end

h = figure;
filename = 'InputTuning';
subplot(2,1,2); hold on;
for i = 1:Ntwk.Input.Source
    plot(xvec, Tuning.Input(i,:)/max(Tuning.Input(:)), 'LineWidth', 2, 'Color', OKeeffe(i,:));
    legendLabels{i} = sprintf('Input %d', i);
end
legend(legendLabels, 'Location', 'best');
xlabel('E neurons'' location (\mum)');
ylabel('Coupling (a.u.)');
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

subplot(2,1,1); hold on;
undersample = 1; %48;
plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'k^', 'MarkerSize', 3);
for i = 1:Ntwk.Input.Source
    Cnncted = any(Ntwk.Cnnct_Input(:,Ntwk.Input.Origins==i),2);
    plot(Ntwk.Exct.Location(Cnncted, 1),...
        Ntwk.Exct.Location(Cnncted, 2), '.', 'Color', OKeeffe(i,:), 'MarkerSize', 6);
    % plot(Ntwk.Input.Location(Ntwk.Input.Origins==i,1), Ntwk.Input.Location(Ntwk.Input.Origins==i,2), '.', 'Color', OKeeffe(i,:), 'MarkerSize', 3);
end
%xlabel('\mum');
ylabel('\mum');
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