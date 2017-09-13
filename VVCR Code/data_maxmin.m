function [Ret] = data_maxmin (Dat);
% Use findpeaks() to determine extrema of dP/dt.

% The peaks must exceed 100 [mmHg/s] and be
% separated by the number of total seconds of the sample*2-1.
[Ret.dPmaxVal, Ret.dPmaxIdx] = findpeaks(Dat.dPdt, 'MinPeakHeight',100, ...
    'MinPeakDistance',length(Dat.dPdt)/(Dat.time_end*2+6));

% find peaks of 'flipped' (negated) data
DataInv = (-1)*Dat.dPdt;

[Ret.dPminVal, Ret.dPminIdx] = findpeaks(DataInv, 'MinPeakHeight',100, ...
    'MinPeakDistance',length(Dat.dPdt)/(Dat.time_end*2+6));
Ret.dPminVal = -Ret.dPminVal;

% make sure only to analyze complete waveforms:
% data must begin with a dP/dt max, and end with dP/dt min
if Ret.dPminIdx(1) < Ret.dPmaxIdx(1)
    Ret.dPminIdx(1) = [];
    Ret.dPminVal(1) = [];
end

if Ret.dPminIdx(end) < Ret.dPmaxIdx(end)
    Ret.dPmaxIdx(end) = [];
    Ret.dPmaxVal(end) = [];
end

% This isn't working... Rvals isn't consistent... Donno why, it should be!
%mysz = length(Ret.dPmaxIdx);
%Ret.Period = zeros(mysz,2);

%Dat.Rvals'
%Ret.dPmaxIdx'
% Find cycle periods (even if dPmax values are skipped)
%for i = 1: mysz
%    [~,St] = find (Dat.Rvals < Ret.dPmaxIdx(i))
%    [~,En] = find (Dat.Rvals > Ret.dPmaxIdx(i))
%    Ret.Period(i,1) = St(end);
%    Ret.Period(i,2) = En(1);
%end

% Find cycle periods (even if dPmax values are skipped)
cyclenIdx = mean([median(diff(Ret.dPminIdx)) median(diff(Ret.dPmaxIdx))]);
Ret.time_per = cyclenIdx*Dat.time_step;
disp(['    data_maxmin: Average Period = ' num2str(Ret.time_per, '%05.3f') ...
    ' sec']);
