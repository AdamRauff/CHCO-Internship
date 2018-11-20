function [Ret] = data_filter (Pat, Pres, dPdt, Rvals)
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
Ret.time_step = Pat.tstp;
Ret.dat_typ = Pat.type;

Wp = Ret.time_step*2*n;
[b,a] = butter(n,Wp);

% Store original (WITT versions) of pressure & dP/dt in WITT variables
% (filtered values are returned below in vars that have no other term and
% are assumed to be used unless filtering somehow is a problem).
Ret.WITTPres = Pres;
Ret.WITTdPdt = dPdt;

% Compute better (6th order central difference) dPdt and d2Pt from base
% WITT data. These are the "unfiltered" original variants. NEVER EVER use
% 1st order differences for pressure acceleration.
Ret.OrigPres = Pres;
[Ret.OrigdPdt, Ret.OrigdP2t] = data_centdiff(Pat.tstp, Pres);

% The dPdt we actually (hope to) use is the filtered "Orig" dP/dt from
% above. Then Pres and d2P/dt2 are the integral and derivative of this
% filtered result (to maintain consistency). Pres(1) is added to the result
% of cumtrapz (constant of integration).
Ret.dPdt = filtfilt(b,a,Ret.OrigdPdt);
Ret.Pres = cumtrapz(Ret.dPdt)*Ret.time_step+Pres(1);
Ret.dP2t = data_centdiff(Pat.tstp, Ret.dPdt);

%% construct time array (4ms from catheter, 1ms from calf pressure DAQ)
Ret.time_end = Ret.time_step*size(Pres,1);
Ret.Time = [Ret.time_step:Ret.time_step:Ret.time_end]';