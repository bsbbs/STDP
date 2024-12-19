function evs = Generator(Ntrial)
%% in unit of seconds
dur.ev1 = .5; % durations of events
dur.ev2 = .5;
% - ITI before event 1, uniform distributed
maxITI = .7;
minITI = .1;
rng(2011);
ITI = rand(Ntrial,1)*(maxITI - minITI) + minITI;
ITI(1) = 0.1;
% duration of event 1
ev1 = ones(Ntrial,1)*dur.ev1;

% - ISI between event 1 and event 2, uniformly distributed
maxISI = .7;
minISI = .1;
rng(2012);
ISI = rand(Ntrial,1)*(maxISI - minISI) + minISI;
% - durations of Event 2
ev2 = ones(Ntrial,1)*dur.ev2;
tp1 = cumsum(ITI + [0;  ev1(1:end-1) + ISI(1:end-1) + ev2(1:end-1)]); % time point of event 1
tp2 = cumsum(ITI + ev1 + ISI + [0; ev2(1:end-1)]); % time point of event 2
evs = [tp1, tp2];

% 
% % creat event boxcar function
% Tdur = ceil(tp2(end) + 3);
% tgrid = 0:tprslt:Tdur;
% r1 = zeros(1+Tdur/tprslt,1);
% r1p = zeros(1+Tdur/tprslt,1);
% r2 = zeros(1+Tdur/tprslt,1);
% for i=1:numel(tp1)
%     r1(tgrid >= tp1(i) & tgrid <= tp1(i)+dur.representation) = values(i,1); % create boxcar event1
%     r1p(tgrid >= tp1(i) & tgrid <= tp1(i)+dur.representation) = values(i,2); % create boxcar event1p
%     r2(tgrid >= tp2(i) & tgrid <= tp2(i)+dur.choice) = values(i,1); % create boxcar event2
% end
% Seq = [r1, r1p, r2];
