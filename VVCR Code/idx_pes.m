function [ivIdx, ivVal] = idx_pes (idxsz, datsz, Dat, Ext, ...
    ivIdx, ivVal, badcyc);
% Find timings for the occurance of Pes. Note that unlike the method-specific
% landmarks, these eventually get put into the Data structure (in data_double).

disp('    data_idx_pes: finding end systolic pressure indices');

% Not sure how to use these right now. How does this go wrong!?!?!
% badcyc.PD = [];
% badcyc.PP = [];

%% COMPUTE [Pes] TIMINGS - PRESSURE ACCEL METHOD
ivIdx.PesP = zeros(idxsz,1);
ivVal.PesP = zeros(idxsz,1);

for i = 1:idxsz

    % Compute start index.
    ESi = Ext.dPminIdx(i);
    dP2Zero = Dat.dP2t(ESi) + 0.01;

    % Position of end systole: @(PA)min before (dP/dt)min. Step backwards from
    % (dP/dt)min until we reach this point. This is Vanderpool's original
    % proposed location.
    while Dat.dP2t(ESi) < dP2Zero
        dP2Zero = Dat.dP2t(ESi);
        ESi = ESi - 1;
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
    ESEi = ESEi + 1;

    % Only keep method 2 if it's earlier than method 1.
    if ESEi < ESi
        ESi = ESEi;
    end
    
    % assign iv*.PesP values
    ivVal.PesP(i) = Dat.Pres(ESi);
    ivIdx.PesP(i) = ESi;
    
end