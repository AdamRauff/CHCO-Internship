function [Ret] = data_filter (dat_typ, Pres, dPdt, Rvals)

% Digital filter of Pres & dP/dt. May be necessary for multiharmoic fit,
% not yet sure. Below, n = filter order (higher gives more rapid
% attenuation after the passband); TimeStep is sampling rate of clinical
% data; Wp is the normalized cutoff (35Hz in this case). A Butterworth
% filter has zero ripple before the passband then slowly attenuates; this
% seems to do a fine job for what we need.

% Copy Rvals into stucture for later use.
Ret.Rvals = Rvals;

n = 35;
if dat_typ
    Ret.time_step = 1/250;
else
    Ret.time_step = 1/1000;
end
Wp = Ret.time_step*2*n;
[b,a] = butter(n,Wp);

% dPdt_filt is the filtered dP/dt, Pres_filt is the integral of this (to
% maintain consistency). Pres(1) is added to the result of cumtrapz, it's
% the constant of integration. Original pressure & dP/dt are stored in
% _orig variables, filtered values go into "no underscore" vars...
Ret.dPdt = filtfilt(b,a,dPdt);
Ret.Pres = cumtrapz(Ret.dPdt)*Ret.time_step+Pres(1);

Ret.OrigPres = Pres;
Ret.OrigdPdt = dPdt;

% FOR DEBUG COMPARISONS - NONFILTERED DATA RETURNED
%Ret.dPdt = dPdt;
%Ret.Pres = Pres;

% construct time array (4ms from catheter, 1ms from calf pressure DAQ)
Ret.time_end = Ret.time_step*size(Pres,1);
Ret.Time = [Ret.time_step:Ret.time_step:Ret.time_end]';
