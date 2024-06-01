function [FR, timep] = PSTH(SpikeTrace, dt)
Mvwndw = 50; % ms
timesteps = size(SpikeTrace, 1);
Stepswndw = Mvwndw/dt;
Mvsteps = floor(timesteps/Stepswndw);
FR = nan(Mvsteps, size(SpikeTrace,2));
timep = nan(Mvsteps, 1);
i = 0;
for ti = floor(Stepswndw/2):Stepswndw:timesteps
    i = i + 1;
    interval = (ti + 1 - Stepswndw/2):(ti + Stepswndw/2);
    FR(i,:) = sum(SpikeTrace(interval,:),1)/Mvwndw*1000;
    timep(i) = ti*dt;
end
