function [ivIdx, ivVal, badcyc] = idx_pes (idxsz, datsz, Dat, Ext, ...
    ivIdx, ivVal);
% Find timings for the occurance of Pes. Note that unlike the method-specific
% landmarks, these eventually get put into the Data structure (in data_double).

disp('    data_idx_pes: finding end systolic pressure indices');

badcyc.P = [];

ivIdx.PesP = zeros(idxsz,1);
ivVal.PesP = zeros(idxsz,1);

%% COMPUTE [Pes] TIMINGS - PRESSURE ACCEL METHOD
for i = 1:idxsz
    
    mybad = 0;

    % Compute start index.
    ESi = Ext.dPminIdx(i);
    dP2Zero = Dat.dP2t(ESi) + 0.01;

    % Position of end systole: @(PA)min before (dP/dt)min. Step backwards from
    % (dP/dt)min until we reach this point. This is Vanderpool's original
    % proposed location.
    while Dat.dP2t(ESi) < dP2Zero
        dP2Zero = Dat.dP2t(ESi);
        ESi = ESi - 1;
        if ESi == 0 && i == 1
            % Ran off front of curve; reset ESi and try again with experimental
            % method immediately below.
            ESi = Ext.dPminIdx(1);
            mybad = 1;
            break;
        end
    end

    ESi = ESi + 1;
    
    % Experimental position of end systole: @(PA)max before (dP/dt)min. Use this
    % from now on (it works!!)
    %
    % Method 1: Step backwards from (PA)min (computed above) until we reach
    % this point. 
    dP2Zero = Dat.dP2t(ESi) - 0.01;
    while Dat.dP2t(ESi) > dP2Zero
        dP2Zero = Dat.dP2t(ESi);
        ESi = ESi - 1;
        if ESi == 0 && i == 1
            % Ran off front of curve again. Let ESi stay at zero (it will be
            % incremented once below, set bad cycle if it still doesn't work.
            mybad = mybad + 1;
            break;
        end
    end

    
    ESi = ESi + 1;

    % Method 2: find absolute minimum of (PA) after RVSP and before (dP/dt)min,
    % using min(). Then step forward from this point. This avoids finding an
    % early (false) (PA)min.
    [~,Pmi] = max(Dat.Pres(Ext.dPmaxIdx(i):Ext.dPminIdx(i)));
    Pmi = Ext.dPmaxIdx(i)+Pmi-1;
    [~,ESEi] = min(Dat.dP2t(Pmi:Ext.dPminIdx(i)));
    ESEi = Pmi+ESEi-1;

    dP2Zero = Dat.dP2t(ESEi) - 0.01;
    while Dat.dP2t(ESEi) > dP2Zero
        dP2Zero = Dat.dP2t(ESEi);
        ESEi = ESEi - 1;
        if ESEi < Pmi
            break;
        end
    end
    
    if ESEi < Pmi && mybad > 0
        disp(['        curve # ' num2str(i, '%02i') ' can''t verify location of Pes.']);
        badcyc.P = [badcyc.P, i];
    end
%   if isoidx_check_bad (i, badcyc.P)
%       continue;
%   end

    ESEi = ESEi + 1;

    % Only keep method 2 if it's earlier than method 1.
    if ESEi < ESi && ESi > 1
        ESi = ESEi;
    elseif ESi == 1
        ESi = ESEi;
    end
    
    % assign iv*.PesP values
    ivVal.PesP(i) = Dat.Pres(ESi);
    ivIdx.PesP(i) = ESi;
    
end