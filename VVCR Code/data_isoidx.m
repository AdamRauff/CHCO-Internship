function [ivIdx, ivVal, badcyc] = data_isoidx (Dat, Ext)
% This code finds indicies and values for the isovolumic sections of the RV
% pressure waveform. There are several approaches for this: Takeuchi, Kind,
% and now, Vanderpool. These are defined as:
% Takeuchi:   IsoC: 20% of (dP/dt)max to (dP/dt)max
%             IsoR: (dP/dt)min to point at which pressure recovers to value
%                   at IsoC start.
% Kind:       IsoC: 20% of (dP/dt)max to 110% of pressure @ (dP/dt)max
%             IsoR: ±20% of pressure @ (dP/dt)min
% Vanderpool: IsoC: Max of (d2P/dt2) just prior to (dP/dt)max to (dP/dt)max
%             IsoR: Same as Takeuchi.
%
% Additionally, we find the end-systolic point, to compute Pes. This can be
% done two ways:
% Hack:       32ms before (dP/dt)min (from dog paper!)
% Vanderpool: Min of (d2P/dt2) just prior to (dP/dt)min.

%% Find isovolumic timings for Takeuchi method points.
disp('    data_isoidx: finding Takeuchi indices');

% Pre allocate variables. Variable name codes are:
%   ivVal - values;
%   ivIdx - indices;
%   1  - 1st method (Takeuchi);
%   2  - 2nd method (Kind);
%   3  - 3rd method (Vanderpool);
%   Ps - positive iso start;
%   Pe - positive iso end;
%   Ns - negative iso start;
%   Ne - negative iso end;

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
    
    % COMPUTE [Ps] TIMINGS
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
            badcyc.T = [badcyc.T, 1]; % add first to list of bad curves
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

    % COMPUTE [Ne] TIMINGS
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
            badcyc.T = [badcyc.T, -i]; % add to list of bad curves
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

               badcyc.T = [badcyc.T, -i]; % add to list of bad curves
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
            badcyc.T = [badcyc.T, -i];
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

% I guess originally I just avoided using the first cycle if it was bad,
% but really I must avoid using any Ps1 values that are bad - those are the
% positive badcyc.T values.
[~,idx] = find(badcyc.T > 0);
if isempty(idx)
    badcyc.K = [];
else
%   badcyc.K = 1;
    badcyc.K = badcyc.T(idx);
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

%% Find Vanderpool times based on pressure acceleration extrema
% (PA)max before (dP/dt)max, end diastole; and (PA)min before (dP/dt)min,
% end systole. The first of these gives us the timing [Ps], while the 
% second allows a non-dog-based calculation of P_es. Like the Takeuchi
% method, the timings of [Pe] and [Ns] come directly from the dP/dt
% extrema, so creating copies of the dPmXX values is all that is needed.

disp('    data_isoidx: finding Vanderpool indices');

ivIdx.Ps3 = zeros(mysz,1);
ivVal.Ps3 = zeros(mysz,1); 
ivIdx.Ne3 = ivIdx.Ne1;
ivVal.Ne3 = ivVal.Ne1;

ivIdx.dPmax3 = Ext.dPmaxIdx;
ivVal.dPmax3 = Ext.dPmaxVal;
ivIdx.dPmin3 = Ext.dPminIdx;
ivVal.dPmin3 = Ext.dPminVal;

% Here I must avoid using any negative badcyc.T values, which indicate that
% the Takeuchi indices failed on Ne1. Additionally I use the Ps1 value to
% determine a bound for Ps3 (the max of PA should be within 2*Ps1 or so);
% so if Ps1 was not found (value of badcyc.T > 0)... well we can still
% attempt with this one, we'll just be more careful!
[~,idx] = find(badcyc.T < 0);
if isempty(idx)
    badcyc.V = [];
else
    badcyc.V = badcyc.T(idx);
end

step = 2*(ivIdx.dPmax1-ivIdx.Ps1);

for i = 1:mysz
    % Bound ED between dPmax and either double the distance between the
    % determined Ps1 value, or more broadly the beginning of this cycle.
    
    % Find range of indices in which to search for (PA)max
    Pend = ivIdx.dPmax3(i);
    
    idx = find(badcyc.T==i);
    if ~isempty(idx) | step(i) <= 0
        if i == 1
            % ~fifth of a period if we're on the first cycle
            step(i) = round(0.2*Dat.time_per/Dat.time_step);
        else
            % ~third of the way between dPmax and previous dPmin. Note that
            % if there are excuded beats this will go very wrong...!
            step(i) = round((ivIdx.dPmax3(i)-ivIdx.dPmin3(i-1))/3);
        end
        disp(['BOUNDING REGION FOR VANDERPOOL Ps INDEX IS NEGATIVE FOR ' ...
            'CYCLE ' num2str(i)]);
        disp(['REGION INDICIES APPROXIMATED AS ' num2str(step(i)) '; ' ...
            'PLEASE REPORT THIS PRESSURE FILE TO DEVELOPER']);
        disp(['IF REGION IS >> ' num2str(round(0.2*Dat.time_per/ ...
            Dat.time_step)) ' USE RESULTS WITH CAUTION']);
    end
    Pstr = ivIdx.dPmax3(i)-step(i);
    if Pstr < 1
        Pstr = 1;
    end
    
    [~, idx] = max(Dat.dP2t(Pstr:Pend));
    ivVal.Ps3(i) = Dat.Pres(Pstr+idx);
    ivIdx.Ps3(i) = Pstr+idx;
    
    % Here, main concern would be that we found the max at the beginning of
    % the interval, which would indicate that the interval needs to be
    % larger (or that, for the case of Pstr = 1, we don't have the max in
    % the dataset).  If idx = 1, then we reject this set.
    if idx == 1
        badcyc.V = [badcyc.V, i];
    end
     
    % We use negative iso values from Takeuchi, so we don't need to find
    % them again. However, this technique has a more robust way of finding
    % Pes, and this is the place to do it. Here, we bound ES to be between
    % dPmin (end) and half the distance between dPmax & dPmin (start).
    Pstr = round(0.5*(ivIdx.dPmin3(i)+ivIdx.dPmax3(i)));
    Pend = ivIdx.dPmin3(i);
        
    % The largest negative value of PA is then the time of end systole.
    [~, idx] = min(Dat.dP2t(Pstr:Pend));
    ivVal.Pes3(i) = Dat.Pres(Pstr+idx);
    ivIdx.Pes3(i) = Pstr+idx;
    
    % Note that if we have a bound region between dPmax & dPmin then the
    % above should be (nearly) foolproof. No checks needed... HOPEFULLY
    
end

%% Remove bad cycles from ivVal, ivIdx vectors.
badcyc.T = sort(unique(abs(badcyc.T)));
badcyc.K = sort(unique(badcyc.K));
badcyc.V = sort(unique(badcyc.V));

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

if ~isempty(badcyc.V)
    for i = length(badcyc.T):-1:1
        j = badcyc.T(i);
        ivVal.Ps3(j) = []; ivIdx.Ps3(j) = [];
        ivVal.Ne3(j) = []; ivIdx.Ne3(j) = [];

        ivVal.dPmax3(j) = []; ivIdx.dPmax3(j) = [];
        ivVal.dPmin3(j) = []; ivIdx.dPmin3(j) = [];
    end
end


disp(['    data_isoidx: ' num2str(mysz,'%02i') ' extrema sets, ' ...
    num2str(length(ivVal.Ps1),'%02i') ' Takeuchi cycles, ' ...
    num2str(length(ivVal.Ps2),'%02i') ' Kind cycles, ', ...
    num2str(length(ivVal.Ps3),'%02i') ' Vanderpool cycles']);
