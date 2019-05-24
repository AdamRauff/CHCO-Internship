function [ivIdx, ivVal, badcyc] = test_isoidx_takeuchi (idxsz, datsz, Dat, ...
    Ext, ivIdx, ivVal)
% Find isovolumic timings for Takeuchi method points.

disp('    data_isoidx_t: finding Takeuchi indices');

ivIdx.Ps1 = zeros(idxsz,1);
ivVal.Ps1 = zeros(idxsz,1);
ivIdx.Ne1 = zeros(idxsz,1);
ivVal.Ne1 = zeros(idxsz,1);

ivIdx.dPmax1 = Ext.dPmaxIdx;
ivVal.dPmax1 = Ext.dPmaxVal;
ivIdx.dPmin1 = Ext.dPminIdx;
ivVal.dPmin1 = Ext.dPminVal;

badcyc.T = [];

% scroll through all maxima
for i = 1:idxsz

    %% COMPUTE [Ps] TIMINGS
    EDi = ivIdx.dPmax1(i);
    PresCut = 0.20*ivVal.dPmax1(i);

    % Start of positive isovolumic time (ivVal.Ps1): 20% (dP/dt)max
    % Step backwards from (dP/dt)max until we reach this point.
    while Dat.dPdt(EDi) > PresCut
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % If the first dP/dt max is too early in the data, the pressure wave
            % does not contain enough information to include. So we remove the
            % first maximum & min. Escape the loop by setting EDi to be index of
            % a minimum.

            disp(['        curve # 01, start of isovolumic contraction ' ...
	        'not captured at start of sample, skipping.']);

            EDi = ivIdx.dPmin1(1); 
            badcyc.T = [badcyc.T, 1]; % add first to list of bad curves
            break;
        end
    end

    % sometimes EDi is 1 or 2 time points away from dP/dt max. This is usually
    % due to a step-like shape of the pressure curve, because of the error 
    % associated with the physcial system of the catheter.
    % if EDi is less then or equal to 3 points away from dP/dt max
    if abs(EDi-ivIdx.dPmax1(i)) <= 3 
        EDi = EDi - 1; % bump EDi one time point back

        % Continuation mark is a flag that keeps track of the 3 point behind EDi
        CONT_MARK = true;

        % try while loop again, additng the condition that it must be more 
        % than 4 points away, and the 3 points before EDi must also be below
        % the peak 
        while CONT_MARK
            if Dat.dPdt(EDi) <= 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 1) < 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 2) < 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 3) < 0.20*ivVal.dPmax1(i) && ...
                abs(EDi-ivIdx.dPmax1(i)) > 3

                CONT_MARK = false;
            end
            if Dat.dPdt(EDi) <= 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 1) < 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 2) < 0.20*ivVal.dPmax1(i) && ...
                Dat.dPdt(EDi - 3) < 0.20*ivVal.dPmax1(i) && ...
                abs(EDi-ivIdx.dPmax1(i)) <= 3

                disp(['        curve # ' num2str(i, '%02i') ...
                    ', too few points in isovolumic contration, ' ...
		            'skipping.']);

                badcyc.T = [badcyc.T, i]; % add to list of bad curves
                break;
            end
            EDi = EDi - 1;
        end
    end
    EDi = EDi - 2;
    % Consistent Start time, despite shifting start of pos iso phase.
    ivIdx.Tst(i) = EDi;

    % assign iv*.Ps1 values
    ivVal.Ps1(i) = Dat.Pres(EDi);
    ivIdx.Ps1(i) = EDi;

    %% COMPUTE [Ne] TIMINGS
    % find iv*.Ne1 point on the other side of the pressure wave; this point has
    % the same pressure (=ivVal.Ps1); ivVal.Ne1 - ivVal.Ps1 is negative
    ESi = ivIdx.Ps1(i)+15;

    % round values to nearest tenth
    while round(Dat.Pres(ESi),1) > round(ivVal.Ps1(i),1)
        ESi = ESi+1;

        if ESi == datsz
            % the last Dat.dPdt min in the data is part of a pressure wave that
            % is not fully contained in the sample; to correct for this, set ESi
            % to be 10 data points (40 ms) prior to the ivVal.Ps1 to force exit
            % of the while loop.
            disp(['        curve # ' num2str(i, '%02i') ...
	        ', (dP/dt)min not captured at end of sample, skipping.']);

            ESi = ivIdx.Ps1(i)-10;
            badcyc.T = [badcyc.T, -i]; % add to list of bad curves
        end
        
        % if algorithm unable to find the negative ivVal.Ps1, it continues to
        % search into proceeding waveforms, cut it off and get rid of that
        % waveform
        
        % first attain sign of dP/dt
        TmpSign = sign(Dat.dPdt(ESi));
        
        % if the derivative is positive, and the time is past dP/dt min (+10
        % time increments), erase waveform
        if TmpSign == 1 
           if ESi-ivIdx.dPmin1(i) > 10
               % set ESi to be 10 data points (40 ms) prior to the ivVal.Ps1 to
               % force exit while loop

               ESi = ivIdx.Ps1(i)-5;
               disp(['        curve # ' num2str(i, '%02i') ...
	           ', can''t find end isovolumic relaxation pressure, ' ...
		   'skipping.']);

               badcyc.T = [badcyc.T, -i]; % add to list of bad curves
           end
        end
    end

    % find which point is closest to the pressure value of ivVal.Ps1
    tempNeg_ivVal.Ps1s = ESi-3:1:ESi+3; % create local neighborhood
    
    % if the last pressure wave form is being evaluated
    if i == idxsz
        
        % check that the temporary neighborhood does not exceed the size of the
        % pressure / time vector
        while max(tempNeg_ivVal.Ps1s) > datsz
            tempNeg_ivVal.Ps1s(end) = [];
        end
    end
    
    % calculate the difference between the local neighborhood
    % (tempNeg_ivVal.Ps1s) and ivVal.Ps1
    Ps1_Diffs = abs(ivVal.Ps1(i)-Dat.Pres(tempNeg_ivVal.Ps1s));

    % find minimum difference
    [~, tempInds] = min(Ps1_Diffs);

    % assign "Negative" iv*.Ps1 values, i.e. iv*.Ne1 values
    ivIdx.Ne1(i) = tempNeg_ivVal.Ps1s(tempInds);
    ivVal.Ne1(i) = Dat.Pres(tempNeg_ivVal.Ps1s(tempInds));

    % if ivIdx.Ne1 < ivIdx.dPmin, or ivIdx.Ne1 is only about 4 points ahead of
    % the min, then we got a problem: there are basically no isovolumic points
    % on the negative side of the curve. We must reject these examples...
    if ivIdx.Ne1(i) <= ivIdx.dPmin1(i) || ...
        (ivIdx.Ne1(i) > ivIdx.dPmin1(i) && ...
        abs(ivIdx.Ne1(i)-ivIdx.dPmin1(i)) <= 3)

        disp(['        curve # ' num2str(i, '%02i') ...
            ', end diastole leads (dP/dt)min, skipping.']);

        % get rid of curve if it is not already marked
        if isempty(find(badcyc.T==i,1))
            badcyc.T = [badcyc.T, -i];
        end

        % ask Hunter about this scenario - Hunter look into this scenario!!! see
        % patient HA000251.&05 (2nd waveform) for example
    end
    
end
