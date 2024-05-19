function InputTuning(Ntwk, Mycolormap, plotdir)
%% Input tuning
xinterval = Ntwk.Scale/60;
xvec = -Ntwk.Scale:xinterval:Ntwk.Scale;
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
    plot(xvec, Tuning.Input(i,:)/max(Tuning.Input(:)), 'LineWidth', 2, 'Color', Mycolormap(i,:));
    legendLabels{i} = sprintf('Input %d', i);
end
legend(legendLabels, 'Location', 'best');
xlabel('E neurons'' location (\mum)');
ylabel('Coupling (a.u.)');
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);

subplot(2,1,1); hold on;
undersample = 48;
plot(Ntwk.Exct.Location(1:undersample:end,1), Ntwk.Exct.Location(1:undersample:end,2), 'k^', 'MarkerSize', 3);
for i = 1:Ntwk.Input.Source
    plot(Ntwk.Input.Location(Ntwk.Input.Origins==i,1), Ntwk.Input.Location(Ntwk.Input.Origins==i,2), '.', 'Color', Mycolormap(i,:), 'MarkerSize', 3);
end
%xlabel('\mum');
ylabel('\mum');
xlim([-Ntwk.Scale, Ntwk.Scale]);
ylim([-Ntwk.Scale/2, Ntwk.Scale/2]);
mysavefig(h, filename, plotdir, 12, [2.5, 2.5], 2);