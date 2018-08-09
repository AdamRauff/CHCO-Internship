function [ivIdx, ivVal, badcyc] = data_isoidx (mysz, Dat, Ext, ivIdx, ...
    ivVal, badcyc)
% Find Vanderpool times based on pressure acceleration extrema
%
% All four indices must be uniquely found per cycle. Additionally (PA)min
% before (dP/dt)min, allows a "non-dog-based" calculation of P_es.

disp('    data_isoidx: finding Vanderpool indices');

ivIdx.Ps3 = zeros(mysz,1);
ivVal.Ps3 = zeros(mysz,1);
ivVal.Ne3 = zeros(mysz,1);
ivIdx.Ne3 = zeros(mysz,1);

ivIdx.dPmax3 = Ext.dPmaxIdx;
ivVal.dPmax3 = Ext.dPmaxVal;
ivIdx.dPmin3 = Ext.dPminIdx;
ivVal.dPmin3 = Ext.dPminVal;

% Unique initialization for Vanderpool: space to store Pes location.
ivVal.Pes3 = zeros(mysz,1);
ivIdx.Pes3 = zeros(mysz,1);

% These are unique from Takeuchi and Kind so we start with a fresh badcyc.
badcyc.V = [];

% This is an estimate, from Takeuchi cycles, of the isovolumic search
% regions. We can take an approach that doesn't use this, but since the
% Takeuchi method is (fairly) reliable, it's a good first guess. Also
% create a limit variable, one-fifth the size of the average period, that
% we use to reject cycles if the search region is too large.
range = 2*(ivIdx.dPmax1-ivIdx.Ps1);
rnglm = round(0.2*Dat.time_per/Dat.time_step);

for i = 1:mysz
    % Bound ED, or (PA)max-Contraction, between dPmax and either double the
    % distance between the determined Ps1 value, or more broadly the from 
    % the beginning of the cycle. Then bound (PA)min-Contraction using the
    % same size range.
    
    % Custom range flag
    crange = 0;
    
    % Find range of indices in which to search for (PA)max-C
    Pend = ivIdx.dPmax3(i);
    
    idx = find(badcyc.T==i);
    if ~isempty(idx) | range(i) <= 0
        crange = 1;
        if i == 1
            % ~fifth of a period if we're on the first cycle
            range(i) = round(0.2*Dat.time_per/Dat.time_step);
        else
            % ~third of the way between dPmax and previous dPmin. Note that
            % if there are excuded beats this will go very wrong...!
            range(i) = round((ivIdx.dPmax3(i)-ivIdx.dPmin3(i-1))/3);
        end
        
        % Check that custom range isn't too extreme.
        if range(i) > rnglm | range(i) <= 0           
            disp(['        curve # ' num2str(i, '%02i') ', bounding ' ...
                'region for Vanderpool Ps index large or negative, ' ...
                'skipping.']);
            badcyc.V = [badcyc.V, i];
            continue;
        end
    end
    
    Pstr = ivIdx.dPmax3(i)-range(i);
    if Pstr < 1
        Pstr = 1;
    end
    
    % find (PA)max before (dP/dt)max - the time of end diastole.
    if data_isoidx_checkrng (Pstr:Pend, i, 'end diastole')
        continue;
    end
    [~, idx] = max(Dat.dP2t(Pstr:Pend));
    ivVal.Ps3(i) = Dat.Pres(Pstr+idx);
    ivIdx.Ps3(i) = Pstr+idx;
    
    % If we found the max at the beginning of the interval, that (might)
    % indicate that the interval needs to be larger (or that, for the case
    % of Pstr = 1, we don't have the max in the dataset).  If idx = 1, then
    % we reject this set.
    if idx == 1
        disp(['        curve # ' num2str(i, '%02i') ', (PA)max ' ...
	        'found at start of search interval for ED, skipping.']);
        badcyc.V = [badcyc.V, i];
        continue;
    end
     
    % Here, we bound ES to be between dPmin (end) and half the distance
    % between dPmax & dPmin (start).
    Pstr = round(0.5*(ivIdx.dPmin3(i)+ivIdx.dPmax3(i)));
    Pend = ivIdx.dPmin3(i);
        
    % Find (PA)min before (dP/dt)min - the time of end systole.
    if data_isoidx_checkrng (Pstr:Pend, i, 'end systole')
        continue;
    end
    [~, idx] = min(Dat.dP2t(Pstr:Pend));
    ivVal.Pes3(i) = Dat.Pres(Pstr+idx);
    ivIdx.Pes3(i) = Pstr+idx;
    % Note that if we have a bound region between dPmax & dPmin then the
    % above should be (nearly) foolproof.
    
    % Finally, we use the step(i) as in the contraction region to find the
    % end of relaxation.
    Pstr = ivIdx.dPmin3(i);
    if crange
        % Use similar idea for range as with ES bounding, above.
        Pend = round(0.5*(ivIdx.dPmin3(i)-ivIdx.dPmax3(i)));
        
        % Check that custom range isn't too extreme.
        if Pend > rnglm | Pend <= 0
            disp(['        curve # ' num2str(i, '%02i') ', bounding ' ...
                'region for Vanderpool end iso relaxation index large ' ...
                'or negative, skipping.']);
            badcyc.V = [badcyc.V, i];
            continue;
        end
        Pend = Pstr + Pend;
    else

        Pend = Pstr + range(i);
    end
    
    if Pend > length(Dat.dP2t)
        disp(['        curve # ' num2str(i, '%02i') ', (PA)max ' ...
	        'search interval exceeds length of data for end iso ' ...
            'relax, skipping.']);
        badcyc.V = [badcyc.V, i];
        continue;
    end
       
    
    % Find (PA)max after (dP/dt)min - the time of end iso relax.
    if data_isoidx_checkrng (Pstr:Pend, i, 'end iso relaxation')
        continue;
    end
    [~, idx] = max(Dat.dP2t(Pstr:Pend));
    ivVal.Ne3(i) = Dat.Pres(Pstr+idx);
    ivIdx.Ne3(i) = Pstr+idx;
    
    % If we found the max at the end of the interval, that (might)
    % indicate that the interval needs to be larger (or that, for the case
    % of Pend = end, we don't have the max in the dataset).  If idx is this
    % length, then we reject this set.
    if idx == length(Pstr:Pend)
        disp(['        curve # ' num2str(i, '%02i') ', (PA)max ' ...
	        'found at end of search interval for end iso relax, ' ...
            'skipping.']);
        badcyc.V = [badcyc.V, i];
        continue;
    end

end

% END data_isoidx_v
end

%% Auxilliary Function(s)

% --- Test range of input, return error if it's too small.
function [ret] = data_isoidx_checkrng (range, cyc, messag)

ret = 0;
if length(range) < 4
    disp(['        curve # ' num2str(cyc, '%02i') ', bounding ' ...
        'region for ' messag ' too small for good fit, skipping.']);
    ret = 1;
end

end
% --- END data_isoidx_checkrng
