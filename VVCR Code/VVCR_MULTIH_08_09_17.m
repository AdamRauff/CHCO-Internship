function [ Res, Pat ] = VVCR_MULTIH_08_09_17( PathName, FileName)

%% (1) Read in data from the given Pat.FileNam
% determine if Pat.FileNam is from calf or humans to apply apprpriate loadp
% function

% if thrid digit/entry of Pat.FileNamname is numeric == human Pat.FileNam
% otherwise, calf Pat.FileNam --> use calf loadp function
if FileName(1) == 'H' && ischar(FileName(2)) && ~isnan(str2double(FileName(3)));
    [Pres, dPdt, Rvals, Pat.Nam, Pat.MRN, Pat.FileNam, ~, ~] = ...
        loadp_10_10_16(PathName, FileName, 1);
    dat_typ = 1;
else
    [Pres, dPdt, Rvals, Pat.FileNam, ~] = ...
        load_calf_p_5_17_17(PathName, FileName, 100);
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

%% (2) Filter pressure data & create time vector.
[Data_O] = data_filter (dat_typ, Pres, dPdt, Rvals);

%% (3) Determine all indexing for analysis.
for i = 1:2
    % Capture Extrema and record number found.
    [Extr] = data_maxmin (Data_O);
    Data_O.time_per = Extr.time_per;
    OrigTot = (length(Extr.dPminIdx)+length(Extr.dPmaxIdx))/2;

    % Call GUI for assessment of extrema. Build input structure.
    % Pressure and derivative, their extrema
    PeakStr.Data = Data_O;
    PeakStr.Extr = Extr;

    % Images
    Green_Check = imread('check.png');
    PeakStr.Green_Check = Green_Check;

    Red_X = imread('ex.png');
    PeakStr.Red_X = Red_X;

    % call on GUI.
    NoPeaksRet = GUI_No_Peaks_10_10(PeakStr);
    clear PeakStr Green_Check Red_X

    % Interpret return structure.
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

    % Find isovolumic timings for Takaguichi & Kind method.
    [ivIdx, ivVal, badcyc] = data_isoidx (Data_O, Extr);

    % If very few timings were found, filtering may be a problem.
    Found  = length(ivIdx.Ps1)/Res.TotNumWaves;
    Reject = Res.TotNumWaves/OrigTot;
    if (i == 1 && Found < 0.5) || (i == 1 && Reject < 0.5)
        if Found < 0.5
            disp(['VVCR_MULTIH: succesfully gated only ' num2str(Found, ...
                '%4.1f%%') ' of max/min cycles. Retrying with raw data.']);
        else
            disp(['VVCR_MULTIH: user rejected  ' num2str(1-Reject, ...
                '%4.1f%%') ' of max/min cycles. Retrying with raw data.']);
        end
        Data_O.FiltPres = Data_O.Pres;
        Data_O.FiltdPdt = Data_O.dPdt;
        Data_O.Pres = Data_O.OrigPres; 
        Data_O.dPdt = Data_O.OrigdPdt; 
        Data_O = rmfield(Data_O, 'OrigPres');
        Data_O = rmfield(Data_O, 'OrigdPdt');
    else
        break;
    end
        
end

% if there were no good pressure waveforms left, then skip patient
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
ICS.Freq_o = double(((2*pi)*Res.TotNumWaves)/(Data.time_end)); % OLD METHOD
ICS.Freq = 2*pi/Data.time_per;
ICS.Pres = Data.Pres;
ICS.dPmaxIdx = ivIdx.dPmax;
ICS.dPminIdx = ivIdx.dPmin;


[FitT, ivSeg, Plot] = fit_takeuchi (ivSeg, Data, ICS);
[FitO] = fit_takeuchi_o (ivSeg, Data, ICS);
[FitK] = fit_kind (ivSeg, ivIdx, Data, FitT);

% What was average frequency ratio between Takeuchi and Data?
%temp = mean(FitT.RCoef(FitT.BadCyc~=1,:));
%AveFreq = temp(3)/(2*pi);
%[double(((2*pi)*Res.TotNumWaves)/(Data.time_end)) 2*pi/Data.time_per]
%[AveFreq 1/Data.time_per AveFreq*Data.time_per]

% Package all structures for passing to GUI_SINU_FIT
SinuStr.FitT = FitT;
SinuStr.FitO = FitO;
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
    BadCyc  = RetStr.FitT.BadCyc;

    % Initialize return structure to contain ALL fitting data.
    Res.FitT = RetStr.FitT;
    Res.FitO = RetStr.FitO;
    Res.FitK = RetStr.FitK;
    Res.P_es = Data.P_es;
    
    % calculate the number of peaks that were evaluated
    Res.numPeaks = 0;
    for i = 1:length(BadCyc)
        if BadCyc(i) ~= 1
            Res.numPeaks = Res.numPeaks + 1;
        end
    end

    %OKAY! here are the final values and we can FINALLY calculate VVCR.
    PIsoMaxT = RetStr.FitT.PIsoMax;
    PIsoMaxO = RetStr.FitO.PIsoMax;
    PIsoMaxK = RetStr.FitK.RCoef(:,1);

    GOOD_P_es  = Data.P_es(BadCyc~=1);
    GOOD_PmxT = PIsoMaxT(BadCyc~=1);
    GOOD_PmxO = PIsoMaxT(BadCyc~=1);
    GOOD_PmxK = PIsoMaxK(BadCyc~=1);

    % mean, std Pes for the waves that fit well
    Res.P_es_Mean = mean(GOOD_P_es);
    Res.P_es_StD  = std(GOOD_P_es);

    % mean, std P_max for the waves that fit well
    Res.PmaxT_Mean = mean(GOOD_PmxT);
    Res.PmaxT_StD  = std(GOOD_PmxT);
    Res.PmaxO_Mean = mean(GOOD_PmxO);
    Res.PmaxO_StD  = std(GOOD_PmxO);
    Res.PmaxK_Mean = mean(GOOD_PmxK);
    Res.PmaxK_StD  = std(GOOD_PmxK);

    %from Uyen Truongs VVCR paper
    Res.VVCRiT_Mean = GOOD_P_es./(GOOD_PmxT-GOOD_P_es); 
    Res.VVCRiT_StD  = std(Res.VVCRiT_Mean);
    Res.VVCRiT_Mean = mean(Res.VVCRiT_Mean);

    Res.VVCRiO_Mean = GOOD_P_es./(GOOD_PmxO-GOOD_P_es); 
    Res.VVCRiO_StD  = std(Res.VVCRiO_Mean);
    Res.VVCRiO_Mean = mean(Res.VVCRiO_Mean);
    Res.VVCRiO_MeanO = Res.P_es_Mean/(Res.PmaxO_Mean-Res.P_es_Mean);

    Res.VVCRiK_Mean = GOOD_P_es./(GOOD_PmxK-GOOD_P_es); 
    Res.VVCRiK_StD  = std(Res.VVCRiK_Mean);
    Res.VVCRiK_Mean = mean(Res.VVCRiK_Mean);

    % Hunter's 
    Res.VVCRnT_Mean = (GOOD_PmxT./GOOD_P_es)-1;
    Res.VVCRnT_StD  = std(Res.VVCRnT_Mean);
    Res.VVCRnT_Mean = mean(Res.VVCRnT_Mean);

    Res.VVCRnO_Mean = (GOOD_PmxO./GOOD_P_es)-1;
    Res.VVCRnO_StD  = std(Res.VVCRnO_Mean);
    Res.VVCRnO_Mean = mean(Res.VVCRnO_Mean);
    Res.VVCRnO_MeanO = (Res.PmaxO_Mean/Res.P_es_Mean)-1;

    Res.VVCRnK_Mean = (GOOD_PmxK./GOOD_P_es)-1;
    Res.VVCRnK_StD  = std(Res.VVCRnK_Mean);
    Res.VVCRnK_Mean = mean(Res.VVCRnK_Mean);

end

end
