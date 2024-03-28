function Seq = CreateEvents(Ntrial, dt)
%%
tprslt = dt; % temporal resolution
dur.representation = .5;
dur.choice = .5;
dur.feedback = 0;
% - ITI, exponential distributed with mean 4s
meanvalue = 4;
maxITI = 12;
minITI = 1;
pracMeanITI = 3.4965; % practical mean value of exponential distribution when considering truncation
tolerance = .01;
ITI = 0;
Pool = exprnd(pracMeanITI,Ntrial*10,1) + minITI;
Pool = Pool(Pool < maxITI);
while mean(ITI) < meanvalue - tolerance || mean(ITI) > meanvalue + tolerance
    ITI = randsample(Pool,Ntrial);
end

% - Event 1 - presentation of value coding for 4 secs
ev1 = ones(Ntrial,1)*dur.representation;

% - ISI, exponential distributed with mean 3s
meanvalue = 3;
maxISI = 8;
minISI = .2;
pracMeanISI = 4.3796;
tolerance = .01;
ISI = 0;
Pool = exprnd(pracMeanISI,Ntrial*10,1) + minISI;
Pool = Pool(Pool < maxISI);
while mean(ISI) < meanvalue - tolerance || mean(ISI) > meanvalue + tolerance
    ISI = randsample(Pool,Ntrial);
end

% - Event 2 - presentation of choice for 2.5 secs and feedback for .5
% sec
ev2 = ones(Ntrial,1)*(dur.choice + dur.feedback);

tpITI = cumsum([0; ITI(1:end-1) + ev1(1:end-1) + ISI(1:end-1) + ev2(1:end-1)]); % time point of ITI
tp1 = cumsum(ITI + [0;  ev1(1:end-1) + ISI(1:end-1) + ev2(1:end-1)]); % time point of value representation
tpISI = cumsum(ITI + ev1 +  [0; ISI(1:end-1) + ev2(1:end-1)]); % time point of ISI
tp2 = cumsum(ITI + ev1 + ISI + [0; ev2(1:end-1)]); % time point of decision

%         h = figure; hold on;
%         plot(tpITI, zeros(size(tp1)),'b.');
%         plot(tp1, zeros(size(tp1)),'r|');
%         plot(tpISI, zeros(size(tp1)),'r.');
%         plot(tp2, zeros(size(tp1)),'b|');


% creat event boxcar function
Tdur = ceil(tp2(end) + 3);
tgrid = 0:tprslt:Tdur;
r1 = zeros(1+Tdur/tprslt,1);
r2 = zeros(1+Tdur/tprslt,1);
for i=1:numel(tp1)
    r1(tgrid >= tp1(i) & tgrid<=tp1(i)+dur.representation)=1; % create boxcar event1
    r2(tgrid >= tp2(i) & tgrid<=tp2(i)+dur.choice)=1; % create boxcar event2
end
Seq = [r1, r2];
end