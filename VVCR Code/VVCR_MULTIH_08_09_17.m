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
    NoPeaksRet = GUI_GateCheck (PeakStr);
    clear PeakStr Green_Check Red_X

    % Interpret return structure.
    if ~isstruct(NoPeaksRet)

        Res = false;
        Pat = false;
        disp('VVCR_MULTIH: GUI_GateCheck closed.');
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
    Found  = double(length(ivIdx.Ps1))/double(Res.TotNumWaves);
    if i == 1 && Found < 0.5
	quest = ['Very few cycles gated compared to total number of cycles. '...
            ' Keep current (filtered) data, load unfiltered data, or ' ...
            'discard patient?'];
        button = questdlg(quest,'Poor Gating','Keep Current', ...
            'Load Unfiltered', 'Discard Patient', 'Load Unfiltered');
        switch button
            case 'Keep Current'
                break;
            case 'Load Unfiltered'
                disp(['VVCR_MULTIH: succesfully gated only ' num2str(100* ...
		    Found,'%5.2f%%') ' of max/min cycles. Retrying ' ...
                    'with raw data.']);
                Data_O.FiltPres = Data_O.Pres;
                Data_O.FiltdPdt = Data_O.dPdt;
                Data_O.Pres = Data_O.OrigPres; 
                Data_O.dPdt = Data_O.OrigdPdt; 
                Data_O = rmfield(Data_O, 'OrigPres');
                Data_O = rmfield(Data_O, 'OrigdPdt');
            case 'Discard Patient'
                Res = true;
                Pat = true;
                disp('VVCR_MULTIH: patient discarded pre-analysis.');
                return
        end
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

% Package all structures for passing to GUI_FitTakeuchi
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
RetStr = GUI_FitTakeuchi (SinuStr);

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
    BadCycT = RetStr.FitT.BadCyc;
    BadCycO = RetStr.FitO.BadCyc;
    BadCycK = RetStr.FitK.BadCyc;

    AnyGood = BadCycT | BadCycO | BadCycK; 

    % Initialize return structure to contain ALL fitting data.
    Res.FitT = RetStr.FitT;
    Res.FitO = RetStr.FitO;
    Res.FitK = RetStr.FitK;
    Res.P_es = Data.P_es;
    
    % calculate the number of peaks that were evaluated
    Res.numPeaksT = sum(~BadCycT);
    Res.numPeaksO = sum(~BadCycO);
    Res.numPeaksK = sum(~BadCycK);

    %OKAY! here are the final values and we can FINALLY calculate VVCR.
    GOOD_PmxT = RetStr.FitT.PIsoMax(BadCycT~=1);
    GOOD_PmxO = RetStr.FitO.PIsoMax(BadCycO~=1);
    GOOD_PmxK = RetStr.FitK.RCoef(BadCycK~=1,1);

    % mean, std Pes and P_max for the waves that fit well
    Res = compute_MeanStd (Res, Data.P_es(AnyGood), 'P_es'); 
    Res = compute_MeanStd (Res, GOOD_PmxT, 'PmaxT');
    Res = compute_MeanStd (Res, GOOD_PmxO, 'PmaxO');
    Res = compute_MeanStd (Res, GOOD_PmxK, 'PmaxK');

    % Compute UT, KSH VVCR (inverse and normal)
    Res = compute_VVCR (Res, Data.P_es(BadCycT~=1), GOOD_PmxT, 'T');
    Res = compute_VVCR (Res, Data.P_es(BadCycO~=1), GOOD_PmxO, 'O');
    Res = compute_VVCR (Res, Data.P_es(BadCycK~=1), GOOD_PmxK, 'K');

    % Were "Vanderpool Points" added?
    Res.VandT = sum(RetStr.FitT.VCyc);
    Res.VandO = sum(RetStr.FitO.VCyc);
end

% END OF VVCR_MULTIH
end

%% Auxilliary Functions to simplify computation of final values.

function [Out] = compute_MeanStd (In, Var, nam)
% Compute Mean & StD of input value, given a name, put into Out structure.

Out = In;
fieldmean = [nam '_Mean'];
fieldstd  = [nam '_StD'];

Out.(fieldmean) = mean(Var);
Out.(fieldstd)  = std(Var);

end

function [Out] = compute_VVCR (In, Pes, Pmx, nam)
% Compute VVCR Mean & StD of pressures, given a name, put into Out structure.

% Normal Ees/Ea or (Pmx-Pes)/Pes
fieldmean = ['VVCRn' nam '_Mean'];
fieldstd  = ['VVCRn' nam '_StD'];

Out.(fieldmean) = (Pmx./Pes)-1;
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

% Inverse Ea/Ees or Pes/(Pmx-Pes)
Out = In;
fieldmean = ['VVCRi' nam '_Mean'];
fieldstd  = ['VVCRi' nam '_StD'];

Out.(fieldmean) = Pes./(Pmx-Pes);
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

end
