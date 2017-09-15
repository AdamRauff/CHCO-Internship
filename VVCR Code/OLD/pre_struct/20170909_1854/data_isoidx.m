function [ivIdx, ivVal, badcyc] = data_isoidx (Dat, Ext)

% Pre allocate variables. Variable name codes are:
%   ivVal - values;
%   ivIdx - indices;
%   m1 - 1st method (Takaguchi);
%   m2 - 2nd method (Kind);
%   Ps - positive iso start;
%   Pe - positive iso end;
%   Ns - negative iso start;
%   Ne - negative iso end;

mysz = length(Ext.dPmaxIdx);
ivIdx.Ps1 = zeros(mysz,1);
ivVal.Ps1 = zeros(mysz,1);
ivIdx.Ne1 = zeros(mysz,1);
ivVal.Ne1 = zeros(mysz,1);

ivIdx.dPmax = Ext.dPmaxIdx;
ivVal.dPmax = Ext.dPmaxVal;
ivIdx.dPmin = Ext.dPminIdx;
ivVal.dPmin = Ext.dPminVal;

badcyc = [];

%% Find isovolumic timings for Takaguichi method points.
disp('    data_isoidx: finding Takaguichi indices');

% scroll through all maxima
for i = 1:mysz

    EDi = ivIdx.dPmax(i);
    
    % Start of positive isovolumic time (ivVal.Ps1): 20% (dP/dt)max
    % Step backwards from (dP/dt)max until we reach this point.
    while Dat.dPdt(EDi) > 0.20*ivVal.dPmax(i)
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % If the first dP/dt max is too early in the data, the pressure
            % wave does not contain enough information to include. So we remove
            % the first maximum & min. Escape the loop by setting EDi to be 
            % index of a minimum.

            EDi = ivIdx.dPmin(1); 
            badcyc = [badcyc, i]; % add to list of bad curves
        end
    end

    % sometimes EDi is 1 or 2 time points away from dP/dt max. This is usually
    % due to a step-like shape of the pressure curve, because of the error 
    % associated with the physcial system of the catheter.
    % if EDi is less then or equal to 3 points away from dP/dt max
    if abs(EDi-ivIdx.dPmax(i)) <= 3 
       EDi = EDi - 1; % bump EDi one time point back

       % Continuation mark is a flag that keeps track of the 3 point behind EDi
       CONT_MARK = true;

       % try while loop again, additng the condition that it must be more 
       % than 4 points away, and the 3 points before EDi must also be below
       % the peak 
       while CONT_MARK
            if Dat.dPdt(EDi) <= 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 1) < 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 2) < 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 3) < 0.20*ivVal.dPmax(i) && ...
                abs(EDi-ivIdx.dPmax(i)) > 3

                CONT_MARK = false;
            end
            if Dat.dPdt(EDi) <= 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 1) < 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 2) < 0.20*ivVal.dPmax(i) && ...
                Dat.dPdt(EDi - 3) < 0.20*ivVal.dPmax(i) && ...
                abs(EDi-ivIdx.dPmax(i)) <= 3

                badcyc = [badcyc, i]; % add to list of bad curves
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

            ESi = ivIdx.Ps1(i)-10;
            badcyc = [badcyc, i]; % add to list of bad curves
        end
        
        % if algorithm unable to find the negative ivVal.Ps1, it continues to
        % search into proceeding waveforms, cut it off and get rid of that
        % waveform
        
        % first attain sign of dP/dt
        TmpSign = sign(Dat.dPdt(ESi));
        
        % if the derivative is positive, and the time is past dP/dt min 
        % (+10 time increments), erase waveform
        if TmpSign == 1 
           if ESi-ivIdx.dPmin(i) > 10
               % set ESi to be 10 data points (40 ms) prior to the ivVal.Ps1
               % to force exit while loop

               ESi = ivIdx.Ps1(i)-5;
               badcyc = [badcyc, i]; % add to list of bad curves
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
    if ivIdx.Ne1(i) <= ivIdx.dPmin(i) || ...
        (ivIdx.Ne1(i) > ivIdx.dPmin(i) && ...
        abs(ivIdx.Ne1(i)-ivIdx.dPmin(i)) <= 3)

        disp(['VVCR_MULTIH: for curve # ',num2str(i), ...
            ', end diastole leads (dP/dt)min, skipping.']);

        % get rid of curve if it is not already marked
        if isempty(find(badcyc==i,1))
            badcyc = [badcyc, i];
        end

        % ask Hunter about this scenario - Hunter look into this scenario!!!
        % see patient HA000251.&05 (2nd waveform) for example
    end
end

%% Find isovolumic timings for Kind method points.
% Pre allocate variables - add in zeros(,1) here to make it work. Variable
% names are: iv2 - IsoVol points for Kind fit; Pn - positive iso end; Ns -
% negative iso start; Ne - negative iso end; Val - value; Idx - index.
%
% Note that Kind method pos iso starts at the same time as Takaguchi, so we
% just use those values (in iv1Ps___ vectors, found above) when needed.

ivVal.Pe2 = zeros(mysz,1);
ivIdx.Pe2 = zeros(mysz,1);
ivVal.Ns2 = zeros(mysz,1);
ivIdx.Ns2 = zeros(mysz,1);
ivVal.Ne2 = zeros(mysz,1);
ivIdx.Ne2 = zeros(mysz,1);

disp('    data_isoidx: finding Kind indices');

for i = 1:mysz
    % Find pres at (dP/dt)max, then find 10% increase, this is end of Kind
    % pos iso segment. 
    Pcut = 1.10*Dat.Pres(ivIdx.dPmax(i));
    Pend = ivIdx.dPmax(i);
    while Dat.Pres(Pend) < Pcut
        Pend = Pend + 1;

        % Don't let it leave this cycle. If we don't see a 10% rise from
        % (dP/dt)max, simply take the maximum pressure as the end.
        if Pend > ivIdx.dPmin(i)
            [~, Pend] = max(Dat.Pres(ivIdx.dPmax(i):ivIdx.dPmin(i)));
            break;
        end
    end
    ivVal.Pe2(i) = Dat.Pres(Pend);
    ivIdx.Pe2(i) = Pend;

    % Find pres at (dP/dt)min, then find ±20% change, these are start and end
    % of Kind neg iso segment.
    Pcut = 1.20*Dat.Pres(ivIdx.dPmin(i));
    Pstr = ivIdx.dPmin(i);
    Plim = Pstr - 5;
    while Dat.Pres(Pstr) < Pcut
        Pstr = Pstr - 1;

        % Don't let this go too far: add this cycle to "bad list" if we back up
        % to (dP/dt)max (which is really far!) or if we go farther than 5 points
        % away (this may need revision).
        if Pstr < ivIdx.dPmax(i) || Pstr == Plim
            badcyc = [badcyc, i];
            break;
        end

    end
    ivVal.Ns2(i) = Dat.Pres(Pstr);
    ivIdx.Ns2(i) = Pstr;

    Pcut = 0.80*Dat.Pres(ivIdx.dPmin(i));
    Pend = ivIdx.dPmin(i);
    Plim = Pend + 5;
    while Dat.Pres(Pend) > Pcut
        Pend = Pend + 1;

        % Don't let this go too far: add this cycle to "bad list" if we exceed
        % the Takaguchi end (i.e. go beyond the isovolumic portion), we go
        % farther than 5 points away (this may need revision), or, if we're on
        % the last cycle, don't go off the end of the data.
        if i == mysz
            badcyc = [badcyc, i];
            break;
        elseif Pend >= ivIdx.Ne1(i) || Pend == Plim
            badcyc = [badcyc, i];
            break;
        end

    end
    ivVal.Ne2(i) = Dat.Pres(Pend);
    ivIdx.Ne2(i) = Pend;

end


%% Remove bad cycles from iv Val, Idx vectors.
% combine  sure badcyc has no redundencies
badcyc = sort(unique(badcyc));

% Remove bad curves
if ~isempty(badcyc)
    % loop through badcyc vector
    for i = length(badcyc):-1:1
        j = badcyc(i);
        ivVal.Ps1(j) = []; ivIdx.Ps1(j) = [];
        ivVal.Ne1(j) = []; ivIdx.Ne1(j) = [];

        ivVal.Pe2(j) = []; ivIdx.Pe2(j) = [];
        ivVal.Ns2(j) = []; ivIdx.Ns2(j) = [];
        ivVal.Ne2(j) = []; ivIdx.Ne2(j) = [];

        ivVal.dPmax(j) = []; ivIdx.dPmax(j) = [];
        ivVal.dPmin(j) = []; ivIdx.dPmin(j) = [];
    end
end
