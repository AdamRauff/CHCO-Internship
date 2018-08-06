function [Ret] = data_filter (dat_typ, Pres, dPdt, Rvals)
% Digital filter of Pres & dP/dt. May be necessary for multiharmoic fit,
% not yet sure. Below, n = filter order (higher gives more rapid
% attenuation after the passband); TimeStep is sampling rate of clinical
% data; Wp is the normalized cutoff (approximately n Hz). A Butterworth
% filter has zero ripple before the passband then slowly attenuates; this
% seems to do a fine job for what we need.

% For future thought: could do a double 6th order difference, or a 6th
% order 2nd derivative (or similar) to get d2P/dt2, then integrate twice to
% get dP/dt and Pres. Certainly smoothing at the level of acceleration
% would yield smoother curves, but maybe too smooth...

%% Preliminaries
% Copy Rvals into stucture for later use.
Ret.Rvals = Rvals;

n = 15; % number is approximate cutoff in Hz.
if dat_typ
    Ret.time_step = 1/250;
else
    Ret.time_step = 1/1000;
end
Wp = Ret.time_step*2*n;
[b,a] = butter(n,Wp);

% Store original (WITT versions) of pressure & dP/dt in _orig variables
% (filtered values are returned below in "no underscore" vars and are
% assumed to be used). PA is taken as first order forward difference in
% this case.
Ret.OrigPres = Pres;
Ret.OrigdPdt = dPdt;

Ret.OrigdP2t = diff(dPdt);
Ret.OrigdP2t(end+1) = Ret.OrigdP2t(end);
Ret.OrigdP2t = Ret.OrigdP2t/Ret.time_step;

%% Do higher order difference, smoothing, in integration -> cleaner pressure
% Compute 6th order diff of dPdt from original Pres vector. This overwrites
% the WITT version, but who cares.
dPdt = data_centdiff(dat_typ, Pres);

% dPdt_filt is the filtered dP/dt, Pres_filt is the integral of this (to
% maintain consistency). Pres(1) is added to the result of cumtrapz, it's
% the constant of integration. Compute second derivative of pressure using
% 6th order diff, as with dP/dt. 
Ret.dPdt = filtfilt(b,a,dPdt);
Ret.Pres = cumtrapz(Ret.dPdt)*Ret.time_step+Pres(1);
Ret.dP2t = data_centdiff(dat_typ, Ret.dPdt);

%% construct time array (4ms from catheter, 1ms from calf pressure DAQ)
Ret.time_end = Ret.time_step*size(Pres,1);
Ret.Time = [Ret.time_step:Ret.time_step:Ret.time_end]';
