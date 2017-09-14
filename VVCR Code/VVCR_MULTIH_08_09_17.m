function [ Res, Pat ] = VVCR_MULTIH_08_09_17( PathName, FileName)

%% (1) Read in data from the given Pat.FileNam
% determine if Pat.FileNam is from calf or humans to apply apprpriate loadp
% function

% if thrid digit/entry of Pat.FileNamname is numeric == human Pat.FileNam
% otherwise, calf Pat.FileNam --> use calf loadp function
if FileName(1) == 'H' && ischar(FileName(2)) && ~isnan(str2double(FileName(3)));
    [Pres, dPdt, Rvals, Pat.Nam, Pat.MRN, Pat.FileNam, ~, ~]=loadp_10_10_16(PathName,FileName,100);
    dat_typ = 1;
else
    [Pres, dPdt, Rvals, Pat.FileNam, ~]=load_calf_p_5_17_17(PathName,FileName,100);
    dat_typ = 0;
end

% check to see if RV array was NOT found by loaddp
if length(Pres) == 1 && Pres == 0

    % set output vars to true --> skip Pat.FileNam, return to runAll
    Res = false;
    Pat = false;
    disp('VVCR_MULTIH: Loadp did not detect an RV column');
    return

end

%% (2,3) Filter pressure data & create time vector; find dP/dt Extrema.
[Data_O] = data_filter (dat_typ, Pres, dPdt, Rvals);
[Extr] = data_maxmin (Data_O);
Data_O.time_per = Extr.time_per;

%% (4) GUI for Manual Deletion of Mins and Maxs

% The GUI allows the user to manually delete inappropriate minima and
% maxima. In order to do that, The following information must be passed
% onto the GUI:

% Pressure and derivative, their extrema
PeakStr.Data = Data_O;
PeakStr.Extr = Extr;

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

    Res = false;
    Pat = false;
    disp('VVCR_MULTIH: GUI_No_Peaks closed.');
    return

end

% if the exit button has been pressed
if NoPeaksRet.TotNumWaves == false

    Res = false;
    Pat = false;
    disp('VVCR_MULTIH: You chose to exit the analysis');
    disp(['    The Pat.FileNam ', FileName, ' was not evaluated!']);
    return

% if the discard patient button has been pressed
elseif NoPeaksRet.TotNumWaves == true

    Res = true;
    Pat = true;
    disp('VVCR_MULTIH: patient discarded pre-analysis.');
    return

% otherwise
else
    
    % update minima and maxima per user filter GUI
    Extr.dPmaxIdx = NoPeaksRet.Extr.dPmaxIdx;
    Extr.dPmaxVal = NoPeaksRet.Extr.dPmaxVal;

    Extr.dPminIdx = NoPeaksRet.Extr.dPminIdx;
    Extr.dPminVal = NoPeaksRet.Extr.dPminVal;

    % obtain number of total waveforms
    Res.TotNumWaves = NoPeaksRet.TotNumWaves;

end

clear NoPeaksRet

%% (5) Find isovolumic timings for Takaguichi & Kind method.
[ivIdx, ivVal, badcyc] = data_isoidx (Data_O, Extr);

% if there were no good pressure waveforms left, then skip patient
% DO SOMETHING ELSE HERE: FILTER AT DIFFERENT FREQUENCY, TRY NO FILTER, ETC.
if isempty(ivIdx.Ps1)

     % set all output variables to true, return to runAll
    Res = true;
    Pat = true;
    return

end

%% (7) Fourier Series interpolation and finding isovolumic segments for fitting
[Data] = data_double (Data_O, ivIdx);

[ivSeg, ivIdx] = data_isoseg (false, Data, ivIdx);

clear Data_O

%% (8) Fit the data

% frequnecy is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
% ICs structure for first pass - enables individual computation of ICs
%ICS.Freq = double(((2*pi)*Res.TotNumWaves)/(Data.time_end)); % OLD METHOD
ICS.Freq = 2*pi/Data.time_per;
ICS.Pres = Data.Pres;
ICS.dPmaxIdx = ivIdx.dPmax;
ICS.dPminIdx = ivIdx.dPmin;


[FitT, ivSeg, Plot] = fit_takeuchi (ivSeg, Data, ICS);
[FitK] = fit_kind (ivSeg, ivIdx, Data, FitT);

% What was average frequency ratio between Takeuchi and Data?
%temp = mean(FitT.RCoef(FitT.BadCyc~=1,:));
%AveFreq = temp(3)/(2*pi);
%[double(((2*pi)*Res.TotNumWaves)/(Data.time_end)) 2*pi/Data.time_per]
%[AveFreq 1/Data.time_per AveFreq*Data.time_per]

% Package all structures for passing to GUI_SINU_FIT
SinuStr.FitT = FitT;
SinuStr.FitK = FitK;
SinuStr.Data = Data;
SinuStr.Plot = Plot;

SinuStr.ivIdx = ivIdx;
SinuStr.ivVal = ivVal;
SinuStr.ivSeg = ivSeg;

%% (9) Visualize / check the fit
% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
RetStr = GUI_SINU_FIT_08_29_17 (SinuStr);

clear SinuStr

%% (10) cleanup, arrange data to return to runAll
if ~isstruct(RetStr)
    % if the exit button has been pressed
    if RetStr == false
    
        % set all output variables to false, return to runAll
        Res = false;
        Pat = false;
        disp('VVCR_MULTIH: You chose to exit the analysis');
        disp(['    The Pat.FileNam ', FileName, ' was not evaluated!']);
        return

    % if the discard patient button has been pressed
    elseif RetStr == true
    
        % set all output variables to true, return to runAll
        Res = true;
        Pat = true;
        disp('VVCR_MULTIH: patient discarded post-analysis.');
        return

    end
% otherwise
else
    % extract Pmax and the list of well fitted curves from the structure
    % that is returned from GUI
    BadCyc  = RetStr.BadCyc;
    PIsoMax = RetStr.PIsoMax;
    
    % calculate the number of peaks that were evaluated
    Res.numPeaks = 0;
    for i = 1:length(BadCyc)
        if BadCyc(i) ~= 1
            Res.numPeaks = Res.numPeaks + 1;
        end
    end

    % Send whole datasets back to runAll (might use them)... this may be
    % unnecessary, delete if so...
    Res.PIsoMax = PIsoMax;
    Res.Pes     = Data.Pes;
    Res.BadCyc  = BadCyc;

    %OKAY! here are the final values and we can FINALLY calculate VVCR.
    GOOD_Pes = Data.P_es(BadCyc~=1);
    GOOD_Pmx = PIsoMax(BadCyc~=1);

    % mean, std Pes for the waves that fit well
    Res.Pes_Mean = mean(GOOD_Pes);
    Res.Pes_StD  = std(GOOD_Pes);

    % mean, std P_max for the waves that fit well
    Res.Pmax_Mean = mean(GOOD_Pmx);
    Res.Pmax_StD  = std(GOOD_Pmx);

    %from Uyen Truongs VVCR paper
    Res.VVCRinv_Mean = GOOD_Pes./(GOOD_Pmx-GOOD_Pes); 
    Res.VVCRinv_StD  = std(Res.VVCRinv_Mean);
    Res.VVCRinv_Mean = mean(Res.VVCRinv_Mean);

    % Hunter's 
    Res.VVCRnorm_Mean = (GOOD_Pmx./GOOD_Pes)-1;
    Res.VVCRnorm_StD  = std(Res.VVCRnorm_Mean);
    Res.VVCRnorm_Mean = mean(Res.VVCRnorm_Mean);

end

end
