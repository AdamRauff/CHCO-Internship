function [Ret] = data_maxmin (Dat)
% Use findpeaks() to determine extrema of dP/dt. Note that this can grab
% multiple peaks per cycle, thus the need to ship the results to GUI_GateCheck.

% The peaks must exceed 100 [mmHg/s] and be separated by the number of total
% seconds of the sample*2-1.
[Ret.dPmaxVal, Ret.dPmaxIdx] = findpeaks(Dat.dPdt, 'MinPeakHeight',100, ...
    'MinPeakDistance',length(Dat.dPdt)/(Dat.time_end*2+6));

% find peaks of 'flipped' (negated) data
DataInv = (-1)*Dat.dPdt;

[Ret.dPminVal, Ret.dPminIdx] = findpeaks(DataInv, 'MinPeakHeight',100, ...
    'MinPeakDistance',length(Dat.dPdt)/(Dat.time_end*2+6));
Ret.dPminVal = -Ret.dPminVal;

% make sure only to analyze complete waveforms: data must begin with a dP/dt
% max, and end with dP/dt min
if Ret.dPminIdx(1) < Ret.dPmaxIdx(1)
    Ret.dPminIdx(1) = [];
    Ret.dPminVal(1) = [];
end

if Ret.dPminIdx(end) < Ret.dPmaxIdx(end)
    Ret.dPmaxIdx(end) = [];
    Ret.dPmaxVal(end) = [];
end

end
