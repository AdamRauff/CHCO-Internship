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

%% (2,3) Filter pressure data & create time vector; find dP/dt Extrema.
[Data_O] = data_filter (dat_typ, Pres, dPdt);

[Extrema] = data_maxmin (Data_O);

%% (4) GUI for Manual Deletion of Mins and Maxs

% The GUI allows the user to manually delete inappropriate minima and
% maxima. In order to do that, The following information must be passed
% onto the GUI:

% Pressure and derivative, their extrema
PeakStr.Data = Data_O;
PeakStr.Ext = Extrema;

% Images
Green_Check = imread('check.png');
PeakStr.Green_Check = Green_Check;

Red_X = imread('ex.png');
PeakStr.Red_X = Red_X;

% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
NoPeaksRet = GUI_No_Peaks_10_10(PeakStr);

clear PeakStr Green_Check Red_X

if ~isstruct(NoPeaksRet)

    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (false);
    disp('VVCR_MULTIH: GUI_No_Peaks closed.');
    return

end

% if the exit button has been pressed
if NoPeaksRet.Ext.dPmaxIdx == false

    % set all output variables to false, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (false);
    disp('VVCR_MULTIH: You chose to exit the analysis');
    disp(['    The file ', FileName, ' was not evaluated!']);
    return

% if the discard patient button has been pressed
elseif NoPeaksRet.TotNumWaves == true

     % set all output variables to true, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (true);
    return

% otherwise
else
    
    % update minima and maxima per user filter GUI
    Extrema.dPmaxIdx = NoPeaksRet.Ext.dPmaxIdx;
    Extrema.dPmaxVal = NoPeaksRet.Ext.dPmaxVal;

    Extrema.dPminIdx = NoPeaksRet.Ext.dPminIdx;
    Extrema.dPminVal = NoPeaksRet.Ext.dPminVal;

    % obtain number of total waveforms
    TotNumWaves = NoPeaksRet.TotNumWaves;

end

clear NoPeaksRet

%% (5) Find isovolumic timings for Takaguichi & Kind method.
[ivIdx, ivVal, badcyc] = data_isoidx (Data_O, Extrema);

% if there were no good pressure waveforms left, then skip patient
% DO SOMETHING ELSE HERE: FILTER AT DIFFERENT FREQUENCY, TRY NO FILTER, ETC.
if isempty(ivIdx.Ps1)

     % set all output variables to true, return to runAll
    [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
        STD_Pes, STD_PMX, TotNumWaves] = deal (true);
    return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iv1PsIdx = ivIdx.Ps1; iv1PsVal = ivVal.Ps1;
iv1NeIdx = ivIdx.Ne1; iv1NeVal = ivVal.Ne1;
iv2PeIdx = ivIdx.Pe2; iv2PeVal = ivVal.Pe2;
iv2NsIdx = ivIdx.Ns2; iv2NsVal = ivVal.Ns2;
iv2NeIdx = ivIdx.Ne2; iv2NeVal = ivVal.Ne2;

dPmaxIdx = ivIdx.dPmax; dPmaxVal = ivVal.dPmax;
dPminIdx = ivIdx.dPmin; dPminVal = ivVal.dPmin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (7) FOURIER SERIES INTERPOLATION
[Data] = data_double (Data_O, ivIdx);
clear Data_O

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pres = Data.Pres;
dPdt = Data.dPdt;
time = Data.Time;
PresDoub = Data.Pres_D;
timeDoub = Data.Time_D;
P_es = Data.P_es;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (7.5) Obtaining Isovolumetric points
% storing isovolumetric data in structure with two fields:
% First field is positive isovolmetric, storing all the points that lie on
% the left side, the positive slope side of the pressure wave, and the
% second field stores the points on the negative slope side.

[ivSeg, ivIdx] = data_isoseg (false, Data, ivIdx);

%% (8) Fit the data

% frequnecy is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
% ICs structure for first pass - enables individual computation of ICs
ICS.Freq = double(((2*pi)*TotNumWaves)/(Data.time_end));
ICS.Pres = Data.Pres;
ICS.dPmaxIdx = ivIdx.dPmax;
ICS.dPminIdx = ivIdx.dPmin;

[FitStr] = isovol_fit (ivSeg, timeDoub, PresDoub, ICS );

% ALL THESE BELOW MIGHT ALL BE SET IN isovol_fit - INDEXES AND VALULES MAY 
% BE UPDATED AS IT FINDS BAD FITS. TWO STRUCTS, ivSeg AND Plot, ARE SET IN
% isovol_fit, SO CAN'T BE UPDATED HERE. ivSeg IS PASSED IN AND UPDATED, WHILE
% Plot IS CREATED WITHIN isovol_fit.
FitStr.Data  = Data;
FitStr.ivIdx = ivIdx;
FitStr.ivVal = ivVal;


FitStr.Plot  = [];    %PLACEHOLDER, REMOVE WHEN Plot IS SET IN isovol_fit
FitStr.ivSeg = ivSeg; %PLACEHOLDER, REMOVE WHEN ivSeg IS UPDATED IN isovol_fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8.5) Set any needed vars that weren't set in isovol_fit.
% Pressure and Time in doubled format:
% These can be set once, don't need to be done on each call to isovol_fit

% Orig time vector with 1/2 points; used for buttondownFcn
FitStr.Time = time; % passed with Data above

% Used as reference for number of peaks
FitStr.iv1PsVal = iv1PsVal; % passed with ivVal above

% Time of EDP, neg EDP occured; used for buttondownFcn
FitStr.iv1PsIdx_D = ivIdx.Ps1_D; % passed with ivIdx above
FitStr.iv1NeIdx_D = ivIdx.Ne1_D; % passed with ivIdx above
% Peak times; used in buttownDownFcn. These are indexs of the old time vector
FitStr.dPmaxIdx = dPmaxIdx; % passed with ivIdx above
FitStr.dPminIdx = dPminIdx; % passed with ivIdx above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (9) Visualize / check the fit
% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
RetStr = GUI_SINU_FIT_08_29_17 (FitStr);

clear FitStr

%% (10) cleanup, arrange data to return to runAll
if ~isstruct(RetStr)
    % if the exit button has been pressed
    if RetStr == false
    
        % set all output variables to false, return to runAll
        [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
            STD_Pes, STD_PMX, TotNumWaves] = deal (false);
        disp('VVCR_MULTIH: GUI_SINU_FIT closed/exit pressed.');
        return

    % if the discard patient button has been pressed
    elseif RetStr == true
    
         % set all output variables to true, return to runAll
        [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, ...
            STD_Pes, STD_PMX, TotNumWaves] = deal (true);
        disp('VVCR_MULTIH: patient discarded.');
        return

    end
% otherwise
else
    % extract Pmax and the list of well fitted curves from the structure
    % that is returned from GUI
    BadCyc  = RetStr.BadCyc;
    PIsoMax = RetStr.PIsoMax;
    
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
