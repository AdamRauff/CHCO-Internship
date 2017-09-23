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

    % set output vars to true --> skip this file, return to runAll, keep going
    [Res, Pat] = deal (true);
    disp('VVCR_MULTIH: Loadp did not detect an RV column');
    return

end

%% (2) Filter pressure data & create time vector.
[Data_O] = data_filter (dat_typ, Pres, dPdt, Rvals);

%% (3) Determine all indexing for analysis.
Res.TotNumWaves = 0;
for i = 1:2
    % Capture Extrema and record number found.
    [Extr] = data_maxmin (Data_O);
    Data_O.time_per = Extr.time_per;

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
    GateStr = GUI_GateCheck (PeakStr);

    [Res, Ret] = interpret_str (GateStr, 'GUI_GateCheck', FileName, Res);
    if Ret
        return;
    end
   
    % update minima and maxima per user filter GUI
    Extr.dPmaxIdx = GateStr.Extr.dPmaxIdx;
    Extr.dPmaxVal = GateStr.Extr.dPmaxVal;
    Extr.dPminIdx = GateStr.Extr.dPminIdx;
    Extr.dPminVal = GateStr.Extr.dPminVal;

    Res.TotNumWaves = GateStr.TotNumWaves;

    % Find isovolumic timings for Takaguichi & Kind method.
    [ivIdx, ivVal, badcyc] = data_isoidx (Data_O, Extr);

    % If very few timings were found, filtering may be a problem.
    mingood = min([length(ivIdx.Ps1) length(ivIdx.Ps2)]);
    Found  = double(mingood)/double(Res.TotNumWaves);
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
                disp(['VVCR_MULTIH: poor fit gating of max/min cycles. ' ...
                    'Retrying with raw data.']);
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
    Res, Pat = deal(true);
    return

end

%% (7) Fourier Series interpolation and finding isovolumic segments for fitting
[Data] = data_double (Data_O, ivIdx);

[ivSeg, ivIdx] = data_isoseg (false, Data, ivIdx);

%% (8) Get initial fits for display/check in GUIs

% frequency is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
% ICs structure for first pass - enables individual computation of ICs
ICS.Freq_o = double(((2*pi)*Res.TotNumWaves)/(Data.time_end)); % OLD METHOD
ICS.Freq = 2*pi/Data.time_per;
ICS.Pres = Data.Pres;
ICS.dPmaxIdx = ivIdx.dPmax1;
ICS.dPminIdx = ivIdx.dPmin1;

[FitT, ivSeg, PlotT] = fit_takeuchi (ivSeg, Data, ICS);
[FitO] = fit_takeuchi_o (ivSeg, Data, ICS);
[FitK, PlotK] = fit_kind (ivSeg, ivIdx, Data, mean(FitT.PIsoMax));

%% (9) Visualize / check the fits
% Structures for unchaning data and for each fit type
GUIDat.ivIdx = ivIdx; GUIDat.ivVal = ivVal; GUIDat.ivSeg = ivSeg;
GUIDat.Data = Data;

% Call the Takeuchi Fit Check GUI
TStr.Plot = PlotT; TStr.FitT = FitT; TStr.FitO = FitO;
RetT = GUI_FitTakeuchi (TStr, GUIDat);
[Res, Ret] = interpret_str (RetT, 'GUI_FitTakeuchi', FileName, Res);
if Ret
    return;
end

% Call the Kind Fit Check GUI
FitK.MeanTP = mean(RetT.FitT.PIsoMax);
KStr.Plot = PlotK; KStr.FitK = FitK;
RetK = GUI_FitKind (KStr, GUIDat);
[Res, Ret] = interpret_str (RetK, 'GUI_FitKind', FileName, Res);
if Ret
    return;
end

RetStr = RetT;
RetStr.FitK = RetK.FitK;

%% (10) cleanup, arrange data to return to runAll
% extract Pmax and the list of well fitted curves from the structure
% that is returned from GUI
BadCycT = RetStr.FitT.BadCyc;
BadCycO = RetStr.FitO.BadCyc;
BadCycK = RetStr.FitK.BadCyc;

% Initialize return structure to contain ALL fitting data.
Res.FitT = RetStr.FitT;
Res.FitO = RetStr.FitO;
Res.FitK = RetStr.FitK;
Res.Pes  = unique([Data.Pes1' Data.Pes2']);
    
% calculate the number of peaks that were evaluated
Res.numPeaksT = sum(~BadCycT);
Res.numPeaksO = sum(~BadCycO);
Res.numPeaksK = sum(~BadCycK);

%OKAY! here are the final values and we can FINALLY calculate VVCR.
GOOD_PmxT = RetStr.FitT.PIsoMax(BadCycT~=1);
GOOD_PmxO = RetStr.FitO.PIsoMax(BadCycO~=1);
GOOD_PmxK = RetStr.FitK.RCoef(BadCycK~=1,1);

% mean, std Pes and P_max for the waves that fit well
Res = compute_MeanStd (Res, Res.Pes, 'Pes'); 
Res = compute_MeanStd (Res, GOOD_PmxT, 'PmaxT');
Res = compute_MeanStd (Res, GOOD_PmxO, 'PmaxO');
Res = compute_MeanStd (Res, GOOD_PmxK, 'PmaxK');

% Compute UT, KSH VVCR (inverse and normal)
Res = compute_VVCR (Res, Data.Pes1(BadCycT~=1), GOOD_PmxT, 'T');
Res = compute_VVCR (Res, Data.Pes1(BadCycO~=1), GOOD_PmxO, 'O');
Res = compute_VVCR (Res, Data.Pes2(BadCycK~=1), GOOD_PmxK, 'K');

% Were "Vanderpool Points" added?
Res.VandT = sum(RetStr.FitT.VCyc);
Res.VandO = sum(RetStr.FitO.VCyc);

% END OF VVCR_MULTIH
end

%% Auxilliary Functions to simplify computation of final values.

function [Res, Ret] = interpret_str (str, guinam, patfile, ResIn)

Res = ResIn;
Ret = 1;

% GUI was closed using OS
if ~isstruct(str)
    Res = false;
    disp(['VVCR_MULTIH: ' guinam ' closed.']);
    return
end

% The exit button has been pressed
if str.Exit == false
    Res = false;
    disp('VVCR_MULTIH: You chose to exit the analysis');
    disp(['    The Pat.FileNam ', patfile, ' was not evaluated.']);

% The discard patient button has been pressed
elseif str.Exit == true
    Res = true;
    disp('VVCR_MULTIH: patient discarded pre-analysis.');

% Everything went well
else
    Ret = 0;

end

end


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

Out = In;

% Normal Ees/Ea or (Pmx-Pes)/Pes
fieldmean = ['VVCRn' nam '_Mean'];
fieldstd  = ['VVCRn' nam '_StD'];

Out.(fieldmean) = (Pmx./Pes)-1;
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

% Inverse Ea/Ees or Pes/(Pmx-Pes)
fieldmean = ['VVCRi' nam '_Mean'];
fieldstd  = ['VVCRi' nam '_StD'];

Out.(fieldmean) = Pes./(Pmx-Pes);
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

end
