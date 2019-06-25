function [ivIdx, ivVal, badcyc] = isoidx_vanderpool (idxsz, datsz, Dat, Ext, ...
    ivIdx, ivVal, badcyc)
% Find isovolumic timings and end systole for Takeuchi method using Vanderpool's
% technique.
%
% Note that Pe and Ns are given by extrema and have already been found.
%
% If the data is unfiltered then DA will be so rough that a "walking" approach
% to find extrema won't work. So, first check if the data provided is the
% filtered, or original. If it's the original data, filter DA heavily...

% This causes the code to move back farther from (dP/dt)min to find Pes (to the
% max of (PA) rather than just to the min of (PA)).
ALT_PES = 1;

disp('    isoidx_vanderpool: finding Vanderpool indices for Takeuchi method');

if ~isfield(Dat, 'OrigdPdt')
    disp(['        using filtered pressure acceleration for landmark ' ...
        'finding.']);
    dP2tOrig = Dat.dP2t;
    Dat.dP2t = Dat.FiltdP2t;
end

ivIdx.Ps3 = zeros(idxsz,1);
ivVal.Ps3 = zeros(idxsz,1);
ivVal.Ne3 = zeros(idxsz,1);
ivIdx.Ne3 = zeros(idxsz,1);

ivIdx.dPmax3 = Ext.dPmaxIdx;
ivVal.dPmax3 = Ext.dPmaxVal;
ivIdx.dPmin3 = Ext.dPminIdx;
ivVal.dPmin3 = Ext.dPminVal;

% These are unique from Takeuchi and Kind so we start with a fresh badcyc.
badcyc.V = [];

for i = 1:idxsz

    %% COMPUTE [Ps] TIMINGS
    EDi = ivIdx.dPmax3(i);
    dP2Zero = Dat.dP2t(EDi) - 0.01;

    % Start of positive isovolumic time (ivVal.Ps3): @(PA)max before (dP/dt)max
    % Step backwards from (dP/dt)max until we reach this point.
    while Dat.dP2t(EDi) > dP2Zero
        dP2Zero = Dat.dP2t(EDi);
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % If the first dP/dt max is too early in the data, the pressure wave
            % does not contain enough information to include. So we remove the
            % first maximum & min. Escape the loop by setting EDi to be index of
            % a minimum.
            disp(['        curve # 01, start of isovolumic contraction ' ...
	        'not captured at start of sample, skipping.']);

            EDi = ivIdx.dPmin3(1); 
            badcyc.V = [badcyc.V, 1]; % add first to list of bad curves
            break;
        end
    end
    if isoidx_check_bad (i, badcyc.V)
        continue;
    end
    EDi = EDi + 1;

    % assign iv*.Ps1 values
    ivVal.Ps3(i) = Dat.Pres(EDi);
    ivIdx.Ps3(i) = EDi;
 

    %% COMPUTE [Ne] TIMINGS
    Eir = ivIdx.dPmin3(i);
    dP2Zero = Dat.dP2t(Eir) - 0.01;

    % Position of end iso relax (ivVal.Ns3): @(PA)max after (dP/dt)min.
    % Step forwards from (dP/dt)min until we reach this point.     
    while Dat.dP2t(Eir) > dP2Zero
        dP2Zero = Dat.dP2t(Eir);
        Eir = Eir + 1;

        if Eir == datsz
            % the last end iso relax is out of the dataset, so we must
            % exclude this set.
            disp(['        curve # ' num2str(i, '%02i') ', end iso '...
	            'relaxation not captured at end of sample, skipping.']);

            Eir = datsz - 10; % protect assignment below from overflow
            badcyc.V = [badcyc.V, i]; % add to list of bad curves
            break;
        end

    end
    if isoidx_check_bad (i, badcyc.V)
        continue;
    end
    Eir = Eir - 1;
    
    % assign iv*.Ne values
    ivVal.Ne3(i) = Dat.Pres(Eir);
    ivIdx.Ne3(i) = Eir;
end

%temp_debug (Dat, ivVal, ivIdx);

end

% --- DEBUG CODE: spit out figures w/Vanderpool indices to visually check them.
function temp_debug (Dat, ivVal, ivIdx)

error('this won''t work with [Pes3] right now, fix it first');

figure('Name', 'Pressure');
plot(Dat.Time, Dat.Pres,'b');
hold on;
h(1) = plot(Dat.Time(ivIdx.Ps1), ivVal.Ps1, 'gs');
h(2) = plot(Dat.Time(ivIdx.Ps3), ivVal.Ps3, 'go');
h(3) = plot(Dat.Time(ivIdx.Ne3), ivVal.Ne3, 'ro');
h(4) = plot(Dat.Time(ivIdx.Pes3), ivVal.Pes3, 'kx');
legend(h, 'Ps1', 'Ps3', 'Ns', 'Pes');

figure('Name', 'd2P/dt2');
plot(Dat.Time, Dat.dP2t,'b');
hold on;
h(1) = plot(Dat.Time(ivIdx.dPmax3), Dat.dP2t(ivIdx.dPmax3),'gx');
h(2) = plot(Dat.Time(ivIdx.Ps3), Dat.dP2t(ivIdx.Ps3),'go');
h(3) = plot(Dat.Time(ivIdx.dPmin3), Dat.dP2t(ivIdx.dPmin3),'rx');
h(4) = plot(Dat.Time(ivIdx.Pes3), Dat.dP2t(ivIdx.Pes3),'r+');
h(5) = plot(Dat.Time(ivIdx.Ne3), Dat.dP2t(ivIdx.Ne3),'ro');
legend(h, 'dPmx', 'Ps', 'dPmn', 'Pes', 'Ne');

end
