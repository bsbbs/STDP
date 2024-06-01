%% Visualize the dynamics of neural activities and synpatic changes

%% visualizing example neurons
if ~exist("WEI", "var")
    load(Rsltfile);
end
%%
dt = 1;
leftt = floor(min(find(Seq(1,:)>0, 1 ), find(Seq(2,:)>0, 1 )*dt)/100)*100;
rightt = ceil((max(find(Seq(1,:)>0, 1 ), find(Seq(2,:)>0, 1 ))*dt+500)/100)*100;
duration = size(Seq,2)*dt; % ms
time = [dt:dt:duration]';
timesteps = numel(time);
%% 
h = figure;
filename = 'Example neurons activity';
subplot(5,3,1);
plot(time, Exmpl.ExctV(:,1), 'k-');
xlim([leftt+30, leftt+200]);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,3,4);
tmp = Exmpl.InhbtV(:,1);
tmp(tmp>=0) = 10;
plot(time, tmp, 'r-');
xlim([leftt+30, leftt+200]);
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,3,7); hold on;
plot(time, squeeze(Exmpl.WEI(:,1,1)), 'k-');
xlim([leftt+30, leftt+200]);
xlabel('Time (ms)');
ylabel('wEI');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,3,10); hold on;
plot(time, squeeze(Exmpl.WIE(:,1,1)), 'r-');
xlim([leftt+30, leftt+200]);
xlabel('Time (ms)');
ylabel('wIE');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,3,13); hold on;
plot(time, squeeze(Exmpl.WEE(:,2,1)), 'k-');
xlim([leftt+30, leftt+200]);
xlabel('Time (ms)');
ylabel('wEE');
mysavefig(h, filename, plotdir, 12, [8,7], 1);

subplot(5,2,2); hold on;
[FR, timep] = PSTH(gather([Exmpl.Espikes(:,1), Exmpl.Ispikes(:,1)]), dt);
plot(timep/60000, FR(:,1), 'k-');
xlim([0 duration/60000]);
xlabel('Time (mins)');
ylabel('Firing rate (Hz)');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,2,4); hold on;
plot(timep/60000, FR(:,2), 'r-');
xlim([0 duration/60000]);
xlabel('Time (mins)');
ylabel('Firing rate (Hz)');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,2,6); hold on;
plot(time/60000, squeeze(Exmpl.WEI(:,1,1)), 'k-');
xlim([0 duration/60000]);
xlabel('Time (mins)');
ylabel('wEI');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,2,8); hold on;
plot(time/60000, squeeze(Exmpl.WIE(:,1,1)), 'r-');
xlim([0 duration/60000]);
xlabel('Time (mins)');
ylabel('wIE');
mysavefig(h, filename, plotdir, 12, [8,7], 1);
subplot(5,2,10); hold on;
plot(time/60000, squeeze(Exmpl.WEE(:,2,1)), 'k-');
xlim([0 duration/60000]);
xlabel('Time (mins)');
ylabel('wEE');
mysavefig(h, filename, plotdir, 12, [8,7], 1);

