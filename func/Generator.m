function [Seq, evs] = Generator(Ntrial, values, dt)
%% Pareto distributions
% Set the random seed for reproducibility
rng(2024);
%%
tprslt = dt/1000; % temporal resolution, in sec
dur.representation = .5;
dur.choice = .5;
dur.feedback = 0;
% - ITI, exponential distributed with mean 1s
meanvalue = 1;
maxITI = 1.9;
minITI = .1;
tolerance = .01;
ITI = rand(Ntrial,1)*(maxITI - minITI) + minITI;
% - Event 1 - presentation of value coding for 4 secs
ev1 = ones(Ntrial,1)*dur.representation;

% - ISI, exponential distributed with mean 1s
meanvalue = 1;
maxISI = 1.9;
minISI = .1;
ISI = rand(Ntrial,1)*(maxISI - minISI) + minISI;

% - Event 2 - presentation of choice for 2.5 secs and feedback for .5
% sec
ev2 = ones(Ntrial,1)*(dur.choice + dur.feedback);

% tpITI = cumsum([0; ITI(1:end-1) + ev1(1:end-1) + ISI(1:end-1) + ev2(1:end-1)]); % time point of ITI
tp1 = cumsum(ITI + [0;  ev1(1:end-1) + ISI(1:end-1) + ev2(1:end-1)]); % time point of value representation
% tpISI = cumsum(ITI + ev1 +  [0; ISI(1:end-1) + ev2(1:end-1)]); % time point of ISI
tp2 = cumsum(ITI + ev1 + ISI + [0; ev2(1:end-1)]); % time point of decision

% h = figure; hold on;
% plot(tpITI, zeros(size(tp1)),'b.');
% plot(tp1, zeros(size(tp1)),'r|');
% plot(tpISI, zeros(size(tp1)),'r.');
% plot(tp2, zeros(size(tp1)),'b|');

% creat event boxcar function
Tdur = ceil(tp2(end) + 3);
tgrid = 0:tprslt:Tdur;
r1 = zeros(1+Tdur/tprslt,1);
r1p = zeros(1+Tdur/tprslt,1);
r2 = zeros(1+Tdur/tprslt,1);
for i=1:numel(tp1)
    r1(tgrid >= tp1(i) & tgrid<=tp1(i)+dur.representation)=values(i,1); % create boxcar event1
    r1p(tgrid >= tp1(i) & tgrid<=tp1(i)+dur.representation)=values(i,2); % create boxcar event1p
    r2(tgrid >= tp2(i) & tgrid<=tp2(i)+dur.choice)=values(i,1); % create boxcar event2
end
Seq = [r1, r1p, r2];
evs = [tp1, tp1, tp2];