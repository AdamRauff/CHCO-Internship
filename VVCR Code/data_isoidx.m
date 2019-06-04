function [ivIdx, ivVal] = data_isoidx (Dat, Ext)
% This code finds indicies and values for the isovolumic sections of the RV
% pressure waveform. There are several approaches for this: Takeuchi, Kind,
% and now, Vanderpool. These are defined as (note d2P/dt2 is abbreivated PA
% or "Pressure Acceleration"):
%
% The points on the curves at which these indices are found are as follows:
%   Pes:        Pes is taken as the time of max PA just prior to (dP/dt)min
%
%   Takeuchi:   IsoC: 20% of (dP/dt)max to (dP/dt)max
%               IsoR: (dP/dt)min to point at which pressure recovers to value
%                     at IsoC start.
%   Kind:       IsoC: 20% of (dP/dt)max to 110% of pressure @ (dP/dt)max
%               IsoR: ±20% of pressure @ (dP/dt)min
%   Vanderpool: IsoC: Max of PA just prior to (dP/dt)max to (dP/dt)max
%               IsoR: (dP/dt)min to max of PA just after (dP/dt)min
%               ...I.E. Vanderpool is TOTALLY DIFFERENT from Takauchi. Also
%               should only be applied to Takeuchi fit model at this time,
%               not to Kind equations.
%
% Variable name key (which primarily appear in the three subfunctions) are as
% follows:
%   ivVal - values;
%   ivIdx - indices;
%   Ps  - positive iso start;
%   Pe  - positive iso end;
%   Ns  - negative iso start;
%   Ne  - negative iso end;
%       1 - 1st method (Takeuchi);
%       2 - 2nd method (Kind);
%       3 - 3rd method (Vanderpool);
%   Pes - location of end systole
%       D - "Dog" experimental method
%       P - Pressure accel method.

datsz = length(Dat.Pres);
idxsz = length(Ext.dPmaxIdx);

%% FUTURE CODE NOTES!!!
% This would be more efficient if just the per-cycle data was passed to each
% function and the loop over the cycles was in this routine. That's too much
% work for now, but it should happen eventually.

%% Call functions to get indices and values.
% Index and Iso-Index routines. idx_pes() finds Pes independent from the three
% landmark-finding codes. Pes is distinct from each method.

[ivIdx, ivVal] = idx_pes (idxsz, datsz, Dat, Ext);

[ivIdx, ivVal, badcyc] = isoidx_takeuchi (idxsz, datsz, Dat, Ext, ivIdx, ...
    ivVal);
[ivIdx, ivVal, badcyc] = isoidx_vanderpool (idxsz, datsz, Dat, Ext, ivIdx, ...
    ivVal, badcyc);
[ivIdx, ivVal, badcyc] = isoidx_kind (idxsz, datsz, Dat, Ext, ivIdx, ...
    ivVal, badcyc);

%% Remove bad cycles from ivVal, ivIdx vectors.
% First sort & unique the indices coming from the iso_idx routines.
badcyc.T = sort(unique(abs(badcyc.T)));
badcyc.V = sort(unique(badcyc.V));
badcyc.K = sort(unique(badcyc.K));

% Remove bad curves by passing fields unique to each type of landmarks to
% removal code. Because Pes doesn't really have bad cycles (that we know of at
% this point), each individual badcyc.{} vector can be used to grab the needed
% values later (no need to delete any Pes values).
[ivVal, ivIdx] = clean_isoidx (badcyc.T, ivVal, ivIdx, ...
    char('Ps1', 'Ne1', 'dPmax1', 'dPmin1'));

[ivVal, ivIdx] = clean_isoidx (badcyc.V, ivVal, ivIdx, ...
        char('Ps3', 'Ne3', 'dPmax3', 'dPmin3')); 
    
[ivVal, ivIdx] = clean_isoidx (badcyc.K, ivVal, ivIdx, ...
    char('Ps2', 'Pe2', 'Ns2', 'Ne2', 'dPmax2', 'dPmin2'));

disp(['    data_isoidx: ' num2str(idxsz,'%02i') ' extrema sets, ' ...
    num2str(length(ivVal.Ps1),'%02i') ' Takeuchi cycles, ' ...
    num2str(length(ivVal.Ps3),'%02i') ' Vanderpool cycles, ', ...
    num2str(length(ivVal.Ps2),'%02i') ' Kind cycles']);
