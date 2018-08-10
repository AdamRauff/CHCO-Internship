function [ivIdx, ivVal, badcyc] = data_isoidx (Dat, Ext)
% This code finds indicies and values for the isovolumic sections of the RV
% pressure waveform. There are several approaches for this: Takeuchi, Kind,
% and now, Vanderpool. These are defined as (note d2P/dt2 is abbreivated PA
% or "Pressure Acceleration":
%
% The points on the curves at which these indices are found are as follows:
% Takeuchi:   IsoC: 20% of (dP/dt)max to (dP/dt)max
%             IsoR: (dP/dt)min to point at which pressure recovers to value
%                   at IsoC start.
% Kind:       IsoC: 20% of (dP/dt)max to 110% of pressure @ (dP/dt)max
%             IsoR: ±20% of pressure @ (dP/dt)min
% Vanderpool: IsoC: Max of PA just prior to (dP/dt)max to (dP/dt)max
%             IsoR: (dP/dt)min to max of PA just after (dP/dt)min
%             Pes : Additionally, Pes is taken as the time of min PA just
%                   prior to (dP/dt)min
%             I.E. Vanderpool is TOTALLY DIFFERENT from Takauchi. Also
%             should only be applied to Takeuchi fit model at this time,
%             not to Kind equations.
%
% Variable name codes (which primarily appear in the three subfunctions)
% are as follows:
%   ivVal - values;
%   ivIdx - indices;
%   Ps  - positive iso start;
%   Pe  - positive iso end;
%   Ns  - negative iso start;
%   Ne  - negative iso end;
%   Pes - location of end systole (Vanderpool only)
%       1  - 1st method (Takeuchi);
%       2  - 2nd method (Kind);
%       3  - 3rd method (Vanderpool);

datsz = length(Dat.Pres);
idxsz = length(Ext.dPmaxIdx);

[ivIdx, ivVal, badcyc] = data_isoidx_t (idxsz, datsz, Dat, Ext);
[ivIdx, ivVal, badcyc] = data_isoidx_v (idxsz, datsz, Dat, Ext, ivIdx, ...
    ivVal, badcyc);
[ivIdx, ivVal, badcyc] = data_isoidx_k (idxsz, datsz, Dat, Ext, ivIdx, ...
    ivVal, badcyc);

%% Remove bad cycles from ivVal, ivIdx vectors.
badcyc.T = sort(unique(abs(badcyc.T)));
badcyc.V = sort(unique(badcyc.V));
badcyc.K = sort(unique(badcyc.K));

% Remove bad curves; unique removal for each fit type.
if ~isempty(badcyc.T)
    for i = length(badcyc.T):-1:1
        j = badcyc.T(i);
        ivVal.Ps1(j) = []; ivIdx.Ps1(j) = [];
        ivVal.Ne1(j) = []; ivIdx.Ne1(j) = [];

        ivVal.dPmax1(j) = []; ivIdx.dPmax1(j) = [];
        ivVal.dPmin1(j) = []; ivIdx.dPmin1(j) = [];
    end
end

if ~isempty(badcyc.V)
    for i = length(badcyc.T):-1:1
        j = badcyc.T(i);
        ivVal.Ps3(j) = []; ivIdx.Ps3(j) = [];
        ivVal.Ne3(j) = []; ivIdx.Ne3(j) = [];

        ivVal.dPmax3(j) = []; ivIdx.dPmax3(j) = [];
        ivVal.dPmin3(j) = []; ivIdx.dPmin3(j) = [];
    end
end

if ~isempty(badcyc.K)
    for i = length(badcyc.K):-1:1
        j = badcyc.K(i);
        ivVal.Ps2(j) = []; ivIdx.Ps2(j) = [];
        ivVal.Pe2(j) = []; ivIdx.Pe2(j) = [];
        ivVal.Ns2(j) = []; ivIdx.Ns2(j) = [];
        ivVal.Ne2(j) = []; ivIdx.Ne2(j) = [];

        ivVal.dPmax2(j) = []; ivIdx.dPmax2(j) = [];
        ivVal.dPmin2(j) = []; ivIdx.dPmin2(j) = [];
    end
end

disp(['    data_isoidx: ' num2str(idxsz,'%02i') ' extrema sets, ' ...
    num2str(length(ivVal.Ps1),'%02i') ' Takeuchi cycles, ' ...
    num2str(length(ivVal.Ps3),'%02i') ' Vanderpool cycles, ', ...
    num2str(length(ivVal.Ps2),'%02i') ' Kind cycles']);


