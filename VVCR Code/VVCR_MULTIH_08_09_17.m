function [ AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, STD_Pes, STD_PMX, TotNumWaves] = VVCR_MULTIH_08_09_17( PathName, FileName)

%% (1) Read in data from the given file
% determine if file is from calf or humans to apply apprpriate loadp
% function

% if thrid digit/entry of filename is numeric == human file
% otherwise, calf file --> use calf loadp function
if FileName(1) == 'H' && ischar(FileName(2)) && ~isnan(str2double(FileName(3)));
    [Pres, dPdt, ~, Pnam, Pmrn, file, ~, ~]=loadp_10_10_16(PathName,FileName,100);
    dat_typ = 1;
else
    [Pres, dPdt, Rvals, file, ~]=load_calf_p_5_17_17(PathName,FileName,100);
    dat_typ = 0;
end

% check to see if RV array was NOT found by loaddp
if length(Pres) == 1 && Pres == 0

    % set all outputs to true --> skip file, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (false);
    disp('VVCR_MULTIH: Loadp did not detect an RV column');
    return

end

%% (2) Filter pressure data & create time vector
% Digital filter of Pres & dP/dt. May be necessary for multiharmoic fit,
% not yet sure. Below, n = filter order (higher gives more rapid
% attenuation after the passband); TimeStep is sampling rate of clinical
% data; Wp is the normalized cutoff (35Hz in this case). A Butterworth
% filter has zero ripple before the passband then slowly attenuates; this
% seems to do a fine job for what we need.
n = 35;
if dat_typ
    time_step = 1/250;
else
    time_step = 1/1000;
end
Wp = time_step*2*n;
[b,a] = butter(n,Wp);

% dPdt_filt is the filtered dP/dt, Pres_filt is the integral of this (to
% maintain consistency). Pres(1) is added to the result of cumtrapz, it's
% the constant of integration. Original pressure & dP/dt are stored in
% _orig variables, filtered values go into "no underscore" vars...
dPdt_filt = filtfilt(b,a,dPdt);
Pres_filt = cumtrapz(dPdt_filt)*time_step+Pres(1);

dPdt_orig = dPdt;
Pres_orig = Pres;

dPdt = dPdt_filt;
Pres = Pres_filt;

% construct time array (units of 4 milliseconds from catheter machine)
% **************** -------------------------
% This time array works for the human data
% what is the time interval from calf data????
% **************** -------------------------
time_end=time_step*size(Pres,1);
time=[time_step:time_step:time_end]';

%% (3) find dP/dt maxima and minima.

% flip data. used for finding dPminVal
DataInv = (-1)*dPdt;

% The peaks must exceed 100 [mmHg/s] and be
% separated by the number of total seconds of the sample*2-1.
[dPmaxVal, dPmaxIdx] = findpeaks(dPdt, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+6)); % 

% find peaks of inverted data
[dPminVal, dPminIdx] = findpeaks(DataInv, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+6));
dPminVal = -dPminVal;

% make sure only to analyze complete waveforms:
% data must begin with a dP/dt max, and end with dP/dt min
if dPminIdx(1) < dPmaxIdx(1)
    dPminIdx(1) = [];
    dPminVal(1) = [];
end

if dPminIdx(end) < dPmaxIdx(end)
    dPmaxIdx(end) = [];
    dPmaxVal(end) = [];
end

%% (4) GUI for Manual Deletion of Mins and Maxs

% The GUI allows the user to manually delete inappropriate minima and
% maxima. In order to do that, The following information must be passed
% onto the GUI:

% Pressure and derivative
PeakStruct.Time_D = time;
PeakStruct.Pres_D = Pres;
PeakStruct.dPdt_D = dPdt;

% dPminVal
PeakStruct.dPminIdx = dPminIdx;
PeakStruct.dPminVal = dPminVal;

% Maxima
PeakStruct.dPmaxIdx = dPmaxIdx;
PeakStruct.dPmaxVal = dPmaxVal;

% pictures
Green_Check = imread('check.png');
PeakStruct.Green_Check = Green_Check;

Red_X = imread('ex.png');
PeakStruct.Red_X = Red_X;

% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
PeakStruct2 = GUI_No_Peaks_10_10(PeakStruct);

clear PeakStruct Green_Check Red_X DataInv

% if the exit button has been pressed
if PeakStruct2.dPmaxIdx == false

    % set all output variables to false, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (false);
    disp('VVCR_MULTIH: You chose to exit the analysis');
    disp(['    The file ', FileName, ' was not evaluated!']);
    return

% if the discard patient button has been pressed
elseif PeakStruct2.TotNumWaves == true

     % set all output variables to true, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (true);
    return

    % otherwise
else
    
    % update minima and maxima per user filter GUI
    dPmaxIdx = PeakStruct2.dPmaxIdx;
    dPmaxVal = PeakStruct2.dPmaxVal;

    dPminIdx = PeakStruct2.dPminIdx;
    dPminVal = PeakStruct2.dPminVal;
    
    % obtain number of total waveforms
    TotNumWaves = PeakStruct2.TotNumWaves;
end


%% (5) Find Iso1StVal: 0.2*dP/dt max
    
% Pre allocate variables - add in zeros(,1) here to make it work
mysz = length(dPmaxIdx);
Iso1StVal = zeros(mysz,1);
Iso1StIdx = zeros(mysz,1);
Iso2StVal = zeros(mysz,1);
Iso2StIdx = zeros(mysz,1);

bad_curve = [];

% scroll through all maxima
for i = 1:length(dPmaxIdx)

    EDi = dPmaxIdx(i);
    
    % while dPdt(EDi) > 0.20*dPdt(dPmaxIdx(i))
    while dPdt(EDi) > 0.20*dPmaxVal(i)
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % the first dP/dt max is too early in the data, and the
            % pressure wave does not contain enough information to
            % include. 

            % Remove the First Maximum and Minimum ****
            EDi = dPminIdx(1); % set EDi to be index of a minimum to force exit of while loop
            bad_curve = [bad_curve, i]; % add to list of bad curves
        end
    end

    % sometimes EDi is 1 or 2 time points away from dP/dt max.
    % This is usually to a step-like shape of the pressure curve, because
    % of the error associated with the physcial system of the catheter.
    % if EDi is less then or equal to 3 points away from dP/dt max
    if abs(EDi-dPmaxIdx(i)) <= 3 
       EDi = EDi - 1; % bump EDi one time point back

       % Continuation mark
       CONT_MARK = true; % this is a flag that keeps track of the 3 point behind EDi
       % try while loop again, additng the condition that it must be more than 4 points away, and the 3 points before EDi must also be below the peak 
       while CONT_MARK
            if dPdt(EDi) <= 0.20*dPmaxVal(i) && dPdt(EDi - 1) < 0.20*dPmaxVal(i) && dPdt(EDi - 2) < 0.20*dPmaxVal(i) && dPdt(EDi - 3) < 0.20*dPmaxVal(i) && abs(EDi-dPmaxIdx(i)) > 3
                CONT_MARK = false;
            end
            if dPdt(EDi) <= 0.20*dPmaxVal(i) && dPdt(EDi - 1) < 0.20*dPmaxVal(i) && dPdt(EDi - 2) < 0.20*dPmaxVal(i) && dPdt(EDi - 3) < 0.20*dPmaxVal(i) && abs(EDi-dPmaxIdx(i)) <= 3
                bad_curve = [bad_curve, i]; % add to list of bad curves
                break
            end
            EDi = EDi - 1;
       end
    end

    % assign Iso1StVal values
    Iso1StVal(i) = Pres(EDi);
    Iso1StIdx(i) = EDi;

    % find point on the other side of the pressure wave, with the
    % same height (pressure) as Iso1StVal. Iso2StVal - Iso1StVal Negative
    ESi = Iso1StIdx(i)+15;

    % round values to nearest tenth
    while round(Pres(ESi),1) > round(Iso1StVal(i),1)
        ESi = ESi+1;

        if ESi == length(Pres)
            % the last dPdt min in the data is part of a pressure
            % wave that is not fully contained in the sample

            ESi = Iso1StIdx(i)-10; % set ESi to be 10 data points (40 ms) prior to the Iso1StVal to force exit while loop
            bad_curve = [bad_curve, i]; % add to list of bad curves
        end
        
        % if algorithm unable to find the negative Iso1StVal, and continues to
        % search into proceeding waveforms, cut it off and get rid of that
        % waveform
        
        % first attain sign of dP/dt
        TmpSign = sign(dPdt(ESi));
        
        % if the derivative is positive, and the time is past dP/dt min (+10 time increments),
        % erase waveform
        if TmpSign == 1 
           if ESi-dPminIdx(i) > 10
               bad_curve = [bad_curve, i]; % add to list of bad curves
               ESi = Iso1StIdx(i)-5; % set ESi to be 10 data points (40 ms) prior to the Iso1StVal to force exit while loop
           end
        end
    end

    % find which point is closest to the pressure value of Iso1StVal
    tempNeg_Iso1StVals = ESi-3:1:ESi+3; % create local neighborhood
    
    % if the last pressure wave form is being evaluated
    if i == length(Iso1StVal)
        
        % check that the temporary neighborhood does not exceed the size of
        % the pressure / time vector
        while max(tempNeg_Iso1StVals)>length(Pres)
            tempNeg_Iso1StVals(end) = [];
        end
    end
    
    % calculate the difference between the local neighborhood
    % (tempNeg_Iso1StVals) and Iso1StVal
    Iso1StVal_Diffs = abs(Iso1StVal(i)-Pres(tempNeg_Iso1StVals));

    % find minimum difference
    [~, tempInds] = min(Iso1StVal_Diffs);

    % assign Negative Iso1StVal values
    Iso2StIdx(i) = tempNeg_Iso1StVals(tempInds); % keep in mind this is an indicie of the time vector
    Iso2StVal(i) = Pres(tempNeg_Iso1StVals(tempInds)); % pressure in mmHg

    % if Iso2StIdx < dPminIdx, or Iso2StIdx is only about 4 points ahead of
    % the min, then we got a problem: there are basically no isovolumic points
    % on the negative side of the curve. We must reject these examples...
    if Iso2StIdx(i) <= dPminIdx(i) || ...
        (Iso2StIdx(i) > dPminIdx(i) && abs(Iso2StIdx(i)-dPminIdx(i)) <= 3)

        disp(['VVCR_MULTIH: for curve # ',num2str(i), ', end diastole leads (dP/dt)min, skipping.']);
        % get rid of curve if it is not already marked
        if isempty(find(bad_curve==i,1))
            bad_curve = [bad_curve, i];
        end

        % ask Hunter about this scenario
        % see patient HA000251.&05 (2nd waveform) for example

    end
end

%% (6) remove bad curves from set
tempNum = length(bad_curve);

% make sure bad_curve has no redundencies
bad_curve = sort(unique(bad_curve));

% Remove bad curves
if ~isempty(bad_curve)
    % loop through bad_curves vector
    for i = tempNum:-1:1
        j = bad_curve(i);
        Iso1StVal(j) = [];
        Iso1StIdx(j) = [];
        Iso2StIdx(j) = [];
        Iso2StVal(j) = [];
        dPmaxIdx(j) = [];
        dPmaxVal(j) = [];
        dPminIdx(j) = [];
        dPminVal(j) = [];
    end
end
clear j

% if there were no good pressure waveforms left, then skip patient
if isempty(Iso1StVal)
    
     % set all output variables to true, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (true);
    return

    % otherwise
end

%% (7) FOURIER SERIES INTERPOLATION
% 
% increase the number of data points to provide a tighter fit of the 
% sinusoid that is fitted to the isovolPresumic points
PresDoub = interpft(Pres,length(Pres)*2);
timeDoub = 0.002:0.002:time(end);

%% Pes = 30 ms prior to dp/dt min %%%%%%%%%%%%%
dtmin_30=time(dPminIdx)-0.03;

% Convert time to integer (can't use double as index)
P_esTimes = timeDoub(uint16(dtmin_30/0.002)); 
P_es = PresDoub(uint16(dtmin_30/0.002));

%% Obtaining Isovolumetric points

% storing isovolumetric data in structure with two fields:
% First field is positive isovolmetric, storing all the points that lie on
% the left side, the positive slope side of the pressure wave, and the
% second field stores the points on the negative slope side.

[DatStr] = data_isovol(Iso1StIdx, Iso2StIdx, time, timeDoub, PresDoub, ...
                       dPmaxIdx, dPminIdx, false);
isovolPres = DatStr.P;
isovolTime = DatStr.T;
Iso1StIdx_Doub = DatStr.i1d;
Iso2StIdx_Doub = DatStr.i2d;

%% (8) Fit the data

% frequnecy is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
% ICs structure for first pass - enables individual computation of ICs
ICS.Freq = double(((2*pi)*TotNumWaves)/(time_end));
ICS.Pres = Pres;
ICS.dPmaxIdx = dPmaxIdx;
ICS.dPminIdx = dPminIdx;

[FitStr] = isovol_fit (isovolPres, isovolTime, timeDoub, PresDoub, ICS );

%% (8.5) Set any needed vars that weren't set in isovol_fit.
% Pressure and Time in doubled format:
% These can be set once, don't need to be done on each call to isovol_fit
FitStr.Time_D = timeDoub;
FitStr.Pres_D = PresDoub;

% Orig time vector with 1/2 points; used for buttondownFcn
FitStr.Time = time;

% Used as reference for number of peaks
FitStr.Iso1StVal = Iso1StVal;

% Time of EDP, neg EDP occured; used for buttondownFcn
FitStr.Iso1StIdx_D = Iso1StIdx_Doub;
FitStr.Iso2StIdx_D = Iso2StIdx_Doub;
% Peak times; used in buttownDownFcn. These are indexs of the old time vector
FitStr.dPmaxIdx = dPmaxIdx;
FitStr.dPminIdx = dPminIdx;

%% (9) Visualize / check the fit
% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
SinuStr = GUI_SINU_FIT_08_29_17 (FitStr);

clear FitStr

%% (10) cleanup, arrange data to return to runAll
% if the exit button has been pressed
if size(SinuStr(1).output,1) == 1 && SinuStr(2).output == false
    
    % set all output variables to false, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (false);
    return

% if the discard patient button has been pressed
elseif size(SinuStr(1).output,1) == 1 && SinuStr(2).output == true
    
     % set all output variables to true, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (true);
    return

% otherwise
else
    % extract Pmax and the list of well fitted curves from the structure
    % that is returned from GUI
    BadCyc  = SinuStr(1).output;
    PIsoMax = SinuStr(2).output;
    
    % calculate the number of peaks that were evaluated
    numPeaks = 0;
    for i = 1:length(BadCyc)
        if BadCyc(i) ~= 1
            numPeaks = numPeaks + 1;
        end
    end
    %OKAY! here are the final values and we can FINALLY calculate VVCR.

    % average Pes for the waves that fit well
    AVG_Pes = mean(P_es(BadCyc~=1)); 
    
    % standard deviation of Pes
    STD_Pes = std(P_es(BadCyc~=1));

    % average P_max for the waves that fit well
    AVG_Pmax = mean(PIsoMax(BadCyc~=1)); 
    
    % standard deviation of Pmax
    STD_PMX = std(PIsoMax(BadCyc~=1));

    %from Uyen Truongs VVCR paper
    VVCR_UT = AVG_Pes/(AVG_Pmax-AVG_Pes); 

    % Hunter's 
    VVCR_KH = (AVG_Pmax/AVG_Pes)-1;
end

end
