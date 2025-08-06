function h = ChannelPlot(timep, FR, Locations)
smpl = 40;
gap = round(numel(Locations)/smpl);
FR = gap*FR./max(FR(:));
[~, I] = sort(Locations);
FR = FR + repmat(I', size(FR,1), 1);
h = figure;
plot(timep, FR(:,1:gap:end), '-');
xlabel('Firing rates (Hz)');
ylabel('Neurons sorted to x location');
end