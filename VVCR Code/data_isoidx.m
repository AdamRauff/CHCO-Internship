function [ivIdx, ivVal, badcyc] = data_isoidx (Dat, Ext)

% Pre allocate variables. Variable name codes are:
%   ivVal - values;
%   ivIdx - indices;
%   1  - 1st method (Takeuchi);
%   2  - 2nd method (Kind);
%   Ps - positive iso start;
%   Pe - positive iso end;
%   Ns - negative iso start;
%   Ne - negative iso end;

%% Find isovolumic timings for Takeuchi method points.
disp('    data_isoidx: finding Takeuchi indices');

mysz = length(Ext.dPmaxIdx);
ivIdx.Ps1 = zeros(mysz,1);
ivVal.Ps1 = zeros(mysz,1);
ivIdx.Ne1 = zeros(mysz,1);
ivVal.Ne1 = zeros(mysz,1);

ivIdx.dPmax1 = Ext.dPmaxIdx;
ivVal.dPmax1 = Ext.dPmaxVal;
ivIdx.dPmin1 = Ext.dPminIdx;
ivVal.dPmin1 = Ext.dPminVal;

badcyc.T = [];

% scroll through all maxima
for i = 1:mysz

    EDi = ivIdx.dPmax1(i);
    
    % Start of positive isovolumic time (ivVal.Ps1): 20% (dP/dt)max
    % Step backwards from (dP/dt)max until we reach this point.
    while Dat.dPdt(EDi) > 0.20*ivVal.dPmax1(i)
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % If the first dP/dt max is too early in the data, the pressure
            % wave does not contain enough information to include. So we remove
            % the first maximum & min. Escape the loop by setting EDi to be 
            % index of a minimum.

            disp(['        curve # 01, start of isovolumic contraction ' ...
	        'not captured at start of sample, skipping.']);

            EDi = ivIdx.dPmin1(1); 
            badcyc.T = [badcyc.T, 1]; % add to list of bad curves
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
                break
            end
            EDi = EDi - 1;
       end
    end

    % assign iv*.Ps1 values
    ivVal.Ps1(i) = Dat.Pres(EDi);
    ivIdx.Ps1(i) = EDi;

    % find iv*.Ne1 point on the other side of the pressure wave; this point
    % has the same pressure (=ivVal.Ps1); ivVal.Ne1 - ivVal.Ps1 is negative
    ESi = ivIdx.Ps1(i)+15;

    % round values to nearest tenth
    while round(Dat.Pres(ESi),1) > round(ivVal.Ps1(i),1)
        ESi = ESi+1;

        if ESi == length(Dat.Pres)
            % the last Dat.dPdt min in the data is part of a pressure wave
            % that is not fully contained in the sample; to correct for this,
            % set ESi to be 10 data points (40 ms) prior to the ivVal.Ps1 to
            % force exit of the while loop.
            disp(['        curve # ' num2str(i, '%02i') ...
	        ', (dP/dt)min not captured at end of sample, skipping.']);

            ESi = ivIdx.Ps1(i)-10;
            badcyc.T = [badcyc.T, i]; % add to list of bad curves
        end
        
        % if algorithm unable to find the negative ivVal.Ps1, it continues to
        % search into proceeding waveforms, cut it off and get rid of that
        % waveform
        
        % first attain sign of dP/dt
        TmpSign = sign(Dat.dPdt(ESi));
        
        % if the derivative is positive, and the time is past dP/dt min 
        % (+10 time increments), erase waveform
        if TmpSign == 1 
           if ESi-ivIdx.dPmin1(i) > 10
               % set ESi to be 10 data points (40 ms) prior to the ivVal.Ps1
               % to force exit while loop

               ESi = ivIdx.Ps1(i)-5;
               disp(['        curve # ' num2str(i, '%02i') ...
	           ', can''t find end isovolumic relaxation pressure, ' ...
		   'skipping.']);

               badcyc.T = [badcyc.T, i]; % add to list of bad curves
           end
        end
    end

    % find which point is closest to the pressure value of ivVal.Ps1
    tempNeg_ivVal.Ps1s = ESi-3:1:ESi+3; % create local neighborhood
    
    % if the last pressure wave form is being evaluated
    if i == length(ivVal.Ps1)
        
        % check that the temporary neighborhood does not exceed the size of
        % the pressure / time vector
        while max(tempNeg_ivVal.Ps1s)>length(Dat.Pres)
            tempNeg_ivVal.Ps1s(end) = [];
        end
    end
    
    % calculate the difference between the local neighborhood
    % (tempNeg_ivVal.Ps1s) and ivVal.Ps1
    ivVal.Ps1_Diffs = abs(ivVal.Ps1(i)-Dat.Pres(tempNeg_ivVal.Ps1s));

    % find minimum difference
    [~, tempInds] = min(ivVal.Ps1_Diffs);

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
            badcyc.T = [badcyc.T, i];
        end

        % ask Hunter about this scenario - Hunter look into this scenario!!!
        % see patient HA000251.&05 (2nd waveform) for example
    end
end

%% Find isovolumic timings for Kind method points.
%
% Kind method states that th pos iso starts at R wave, which isn't available,
% so we approximate by using the same time as Takeuchi. However, when waves
% are deleted we need unique vectors for each method, so we copy the Takeuchi
% (iv*.Ps1) vetors into the Kind vectors. We also create unique copies of the
% extrema vectors for the same reason.

disp('    data_isoidx: finding Kind indices');

ivIdx.Ps2 = ivIdx.Ps1;
ivVal.Ps2 = ivVal.Ps1; 
ivVal.Pe2 = zeros(mysz,1);
ivIdx.Pe2 = zeros(mysz,1);
ivVal.Ns2 = zeros(mysz,1);
ivIdx.Ns2 = zeros(mysz,1);
ivVal.Ne2 = zeros(mysz,1);
ivIdx.Ne2 = zeros(mysz,1);

ivIdx.dPmax2 = Ext.dPmaxIdx;
ivVal.dPmax2 = Ext.dPmaxVal;
ivIdx.dPmin2 = Ext.dPminIdx;
ivVal.dPmin2 = Ext.dPminVal;

[~,idx] = find(badcyc.T==1);
if isempty(idx)
    badcyc.K = [];
else
    badcyc.K = [1];
end

% # of points it can look before it stops... each point is 4ms.
lenlim = 10;

for i = 1:mysz
    % Find pres at (dP/dt)max, then find 10% increase, this is end of Kind
    % pos iso segment. 
    Pcut = 1.10*Dat.Pres(ivIdx.dPmax2(i));
    Pend = ivIdx.dPmax2(i);
    while Dat.Pres(Pend) < Pcut
        Pend = Pend + 1;

        % Don't let it leave this cycle. If we don't see a 10% rise from
        % (dP/dt)max, simply take the maximum pressure as the end.
        if Pend > ivIdx.dPmin2(i)
            [~, Pend] = max(Dat.Pres(ivIdx.dPmax2(i):ivIdx.dPmin2(i)));
            break;
        end
    end
    ivVal.Pe2(i) = Dat.Pres(Pend);
    ivIdx.Pe2(i) = Pend;

    % Find pres at (dP/dt)min, then find ±20% change, these are start and end
    % of Kind neg iso segment.
    Pcut = 1.20*Dat.Pres(ivIdx.dPmin2(i));
    Pstr = ivIdx.dPmin2(i);
    Plim = Pstr - lenlim;
    while Dat.Pres(Pstr) < Pcut
        Pstr = Pstr - 1;

        % Don't let this go too far: add this cycle to "bad list" if we back up
        % to (dP/dt)max (which is really far!) or if we go farther than 5 points
        % away (this may need revision).
        if Pstr < ivIdx.dPmax2(i) || Pstr == Plim
            disp(['        curve # ' num2str(i, '%02i') ', can''t find ' ...
                'start of negiso segment, skipping.']);
            badcyc.K = [badcyc.K, i];
            break;
        end

    end
    ivVal.Ns2(i) = Dat.Pres(Pstr);
    ivIdx.Ns2(i) = Pstr;

    Pcut = 0.80*Dat.Pres(ivIdx.dPmin2(i));
    Pend = ivIdx.dPmin2(i);
    Dend = length(Dat.Pres);
    Plim = Pend + lenlim;
    while Dat.Pres(Pend) > Pcut
        Pend = Pend + 1;

        % Don't let this go too far: add this cycle to "bad list" if we exceed
        % the Takeuchi end (i.e. go beyond the isovolumic portion), we go
        % farther than 5 points away (this may need revision), or, if we're on
        % the last cycle, don't go off the end of the data.
        if Pend > Dend
            disp(['        curve # ' num2str(i, '%02i') ', end of negiso ' ...
	        'segment not captured at end of sample, skipping.']);

            Pend = Dend;
            badcyc.K = [badcyc.K, i];
            break;
        elseif Pend >= ivIdx.Ne1(i) || Pend == Plim
            disp(['        curve # ' num2str(i, '%02i') ', can''t find '...
                'end of negiso segment, skipping.']);

            badcyc.K = [badcyc.K, i];
            break;
        end

    end
    ivVal.Ne2(i) = Dat.Pres(Pend);
    ivIdx.Ne2(i) = Pend;

end

%% Remove bad cycles from ivVal, ivIdx vectors.
badcyc.T = sort(unique(badcyc.T));
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

disp(['    data_isoidx: ' num2str(mysz,'%02i') ' extrema sets, ' ...
    num2str(length(ivVal.Ps1),'%02i') ' Takeuchi cycles, ' ...
    num2str(length(ivVal.Ps2),'%02i') ' Kind cycles']);
