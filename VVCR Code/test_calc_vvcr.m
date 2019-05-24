function [ Res, Pat ] = test_calc_vvcr(PathName, FileName, runNum, isopos)
% This code is designed to run multiple passes, GUI-free, over the Tedford data.
% Thus, the GUI logical, and all possibilities of GUI calls, are removed.
% Further, the Tedford data always gets filtered; the code never asks to load
% the unfiltered data.

%% (1) Read in data from the given FileName
% determine if FileName is from calf or humans to apply apprpriate loadp func if
% third digit/entry of FileName is numeric == human; otherwise, calf.
if FileName(1) == 'H' && ischar(FileName(2)) && ~isnan(str2double(FileName(3)))
    [Pres, dPdt, Rvals, Pat] = load_pres_patient(PathName, FileName);
else
    [Pres, dPdt, Rvals, Pat] = load_pres_calf(PathName, FileName);
end

% check to see if RV array was NOT found by loaddp
if length(Pres) == 1

    % set output vars to true --> skip this file, return to runAll, keep going
    [Res, Pat] = deal (true);
    if Pres == 0
        disp('calc_vvcr: Loadp did not detect an RV column');
    else
        disp('calc_vvcr: Loadp could not open file');
    end
    return

end

%% (2) Filter pressure data, find pressure acceleration, & create time vector.
[Data_O] = data_filter (Pat, Pres, dPdt, Rvals);

%% (3) Determine all indexing for analysis.
Res.A_TotNumWaves = 0;

% %%% Looping removed for Tedford problem. ALWAYS use filter. %%%

% Capture Extrema and record number found, then build input structure for
% GUI, then call.
[Extr] = data_maxmin (Data_O);

%%%%%%%%%%% TEDFORD DATA SPECIFIC CODE %%%%%%%%%%%
if length(Extr.dPmaxVal) == 2 && length(Extr.dPminVal) == 1
    Extr.dPmaxVal(2) = [];
    Extr.dPmaxIdx(2) = [];
end
%%%%%%%%%%% TEDFORD DATA SPECIFIC CODE %%%%%%%%%%%

PeakStr.Data = Data_O;
PeakStr.Extr = Extr;

% Images
Green_Check = imread('check.png');
PeakStr.Green_Check = Green_Check;

Red_X = imread('ex.png');
PeakStr.Red_X = Red_X;

% %%% Original place for GUI call %%%
if length(Extr.dPmaxVal) ~= length(Extr.dPminVal)
    error('Extrema Lengths don''t match, clean up the data!');
else
    GateStr = PeakStr;
    GateStr.Exit = 'good';
    GateStr.TotNumWaves = length(Extr.dPmaxVal);
end

% update minima and maxima per user filter GUI
Extr.dPmaxIdx = GateStr.Extr.dPmaxIdx;
Extr.dPmaxVal = GateStr.Extr.dPmaxVal;
Extr.dPminIdx = GateStr.Extr.dPminIdx;
Extr.dPminVal = GateStr.Extr.dPminVal;

Res.A_TotNumWaves = GateStr.TotNumWaves;
clear('GateStr');

% Find cycle periods, AFTER GUI_GateCheck
if length(Extr.dPmaxIdx) > 1 && length(Extr.dPmaxVal) > 1
    cyclenIdx = mean([median(diff(Extr.dPminIdx)) ...
        median(diff(Extr.dPmaxIdx))]);
    Data_O.time_per = cyclenIdx*Data_O.time_step;
else
    Data_O.time_per = Data_O.Time(end);
end
Res.A_HR = 60/Data_O.time_per;

disp(['    calc_vvcr: Average Period = ' num2str(Data_O.time_per, ...
    '%05.3f') ' sec, timestep ' num2str(Data_O.time_step*1000, '%04.1f') ...
    ' ms']);

%% (4) Find isovolumic timings for Takaguichi & Kind method.
[ivIdx, ivVal] = test_data_isoidx (Data_O, Extr);

% If very few timings were found, filtering may be a problem.
% %%% Removed for Tedford Data %%%
mingood = min([length(ivIdx.Ps1) length(ivIdx.Ps2)]);
Found  = double(mingood)/double(Res.A_TotNumWaves);
% %%% End original loop %%%

RunT = logical(length(ivIdx.Ps1));
RunK = logical(length(ivIdx.Ps2));
RunV = logical(length(ivIdx.Ps3));

% if there were no good pressure waveforms left, then skip patient
if ~RunT && ~RunK && ~RunV

    % set all output variables to true, return to runAll
    [Res, Pat] = deal(true);
    return

end

% Only compute Kind method (for now)...
RunT = 0;
RunV = 0;

%% (5) Fourier Series interpolation and finding isovolumic segments for fitting
% Note that Pes recalc (for PesD) and segmentation occurs in data_double.
[Data] = test_data_double (Data_O, ivIdx);

[ivSeg, ivIdx] = data_isoseg (false, Data, ivIdx);

[ivSeg, ivVal] = test_uniq_isoseg (ivSeg, ivVal, ivIdx, Data, isopos, runNum);

% Get Pes processing done prior to fits: this is for all measured Pes, not just
% the per-cycle (accepted) fits... that might be bad? Not sure. Below, when VVCR
% is computed, only Pes for accepted cycles gets passed.
Res = compute_MeanStd (Res, Data.PesD, 'B_PesD'); 
Res = compute_MeanStd (Res, Data.PesP, 'B_PesP');

% In all the fits below, FHEADx values prior to entry to each RunX section
% define the headers for returned fields. The X_ in front of each name allows
% them to be sorted as desired when returned to runAll. Note that other fields
% already set (like the B_PesP above) also have X_ in front to allow sorting.

%% (6) Perform Takeuchi fit(s), put up check GUI, and compute return quantities
% FitT is "new" (time-normalized) w/new Pes, FitO is same fit (Pmax) with old
% dog Pes.
FHEAD1 = 'E_TakeNorm';
FHEAD2 = 'G_TakeOldM';
if RunT
    % frequency is the conversion to angular frequency 2*pi/T multiplied by the
    % number of waves found over the time period ICs structure for first pass -
    % enables individual computation of ICs
    ICS.Freq_o = double(((2*pi)*Res.A_TotNumWaves)/(Data.time_end)); % OLD METHOD
    ICS.Freq = 2*pi/Data.time_per;
    ICS.Pres = Data.Pres;
    ICS.dPmaxIdx = ivIdx.dPmax1;
    ICS.dPminIdx = ivIdx.dPmin1;

    % Keep non-normalized fit for now, but don't include its output anymore.
    % When time allows, remove the call below [FitO] and its processing in the
    % GUI.
    [FitT, ivSeg, PlotT] = fit_takeuchi (ivSeg, Data, ICS, 1);
    [FitO] = fit_takeuchi (ivSeg, Data, ICS, 0);

    % Call the Takeuchi Fit Check GUI
    TStr.FitT = FitT; TStr.FitO = FitO;
    RetT = TStr;
   
    BadCycT = RetT.FitT.BadCyc;
    % BadCycO = RetT.FitO.BadCyc | RetT.FitT.BadCyc;
    
    % Obviously when FitO actually disappears, this below can be simplified.
    Res.([FHEAD1 '_AllDat']) = RetT.FitT;
    Res.([FHEAD2 '_AllDat']) = RetT.FitT;
    Res.([FHEAD1 '_nFit']) = sum(~BadCycT);
    Res.([FHEAD2 '_nFit']) = sum(~BadCycT);
    Res.([FHEAD1 '_Vcorr']) = sum(RetT.FitT.VCyc);
    Res.([FHEAD2 '_Vcorr']) = sum(RetT.FitO.VCyc);

    GOOD_PmxT = RetT.FitT.PIsoMax(BadCycT~=1);
    GOOD_PmxO = RetT.FitO.PIsoMax(BadCycT~=1);
    Res = compute_MeanStd (Res, GOOD_PmxT, [FHEAD1 '_Pmax']);
    Res = compute_MeanStd (Res, GOOD_PmxO, [FHEAD2 '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycT~=1), GOOD_PmxT, FHEAD1);
    Res = compute_VVCR (Res, Data.PesD(BadCycT~=1), GOOD_PmxO, FHEAD2);
else
    [Res] = create_blank_fields (FHEAD1, Res, false);
    [Res] = create_blank_fields (FHEAD2, Res, false);
end

%% (7) Perform Takeuchi fit w/Vanderpool Landmarks, put up check GUI, and 
% compute return quantities. FitV uses "new" ICs (w/fitting limits and new as
% with above. No "old" (unconstrained) fit here.
FHEAD1 = 'F_Vander';
if RunV
    % frequency is the conversion to angular frequency 2*pi/T multiplied by the
    % number of waves found over the time period ICs structure for first pass -
    % enables individual computation of ICs
    ICS.Freq = 2*pi/Data.time_per;
    ICS.Pres = Data.Pres;
    ICS.dPmaxIdx = ivIdx.dPmax3;
    ICS.dPminIdx = ivIdx.dPmin3;

    [FitV, ivSeg, PlotV] = fit_takeuchi (ivSeg, Data, ICS, 2);

    % Call the Vanderpool Fit Check GUI
    VStr.FitV = FitV;
    RetV = VStr;

    BadCycV = RetV.FitV.BadCyc;

    Res.([FHEAD1 '_AllDat']) = RetV.FitV;
    Res.([FHEAD1 '_nFit']) = sum(~BadCycV);
    Res.([FHEAD1 '_Vcorr']) = sum(RetV.FitV.VCyc);
    
    GOOD_PmxV = RetV.FitV.PIsoMax(BadCycV~=1);
    Res = compute_MeanStd (Res, GOOD_PmxV, [FHEAD1 '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycV~=1), GOOD_PmxV, FHEAD1);
else
    [Res] = create_blank_fields (FHEAD1, Res, false);
end

%% (8) Perform Kind fit, put up check GUI, and compute return quantities
% FitK is unweighted residuals with time normalization. FitN has contraction
% error weighted to be (roughly) the same as relaxtion error.
FHEAD1 = 'C_KindNorm';
FHEAD2 = 'D_KindExpr';
if RunK
    % FitN method can also be eliminated (or converted into experimental method
    % for computing 
    [FitK, PlotK] = test_fit_kind (ivSeg, ivIdx, Data, 0);
    [FitN] = test_fit_kind (ivSeg, ivIdx, Data, 1);

    % Call the Kind Fit Check GUI 
    % Not sure why MeanTP was being sent to fit check. Isn't used once it
    % arrives. Comment out for now, remove later if it's never used...   
    %if RunT
    %    FitK.MeanTP = mean(RetT.FitT.PIsoMax);
    %else
    %    FitK.MeanTP = mean(RetT.FitT.PIsoMax);
    %end
    KStr.FitK = FitK;
    RetK = KStr;

    BadCycK = RetK.FitK.BadCyc;
    BadCycN = FitN.BadCyc | RetK.FitK.BadCyc;

    Res.([FHEAD1 '_AllDat']) = RetK.FitK;
    Res.([FHEAD2 '_AllDat']) = FitN;
    Res.([FHEAD1 '_nFit']) = sum(~BadCycK);
    Res.([FHEAD2 '_nFit']) = sum(~BadCycN);

    GOOD_PmxK = RetK.FitK.RCoef(BadCycK~=1,1);
    GOOD_PmxN = FitN.RCoef(BadCycN~=1,1);
    Res = compute_MeanStd (Res, GOOD_PmxK, [FHEAD1 '_Pmax']);
    Res = compute_MeanStd (Res, GOOD_PmxN, [FHEAD2 '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycK~=1), GOOD_PmxK, FHEAD1);
    Res = compute_VVCR (Res, Data.PesP(BadCycN~=1), GOOD_PmxN, FHEAD2);

else
    [Res] = create_blank_fields (FHEAD1, Res, false);
    [Res] = create_blank_fields (FHEAD2, Res, false);
end

Res = compute_MeanStd (Res, ivVal.dPmax2, 'dPmax');
Res = compute_MeanStd (Res, ivVal.dPmin2, 'dPmin');
Res = compute_MeanStd (Res, ivVal.Ps2, 'iso1Ps');
Res = compute_MeanStd (Res, ivVal.Pe2, 'iso2Pe');
Res = compute_MeanStd (Res, ivVal.Ns2, 'iso3Ns');
Res = compute_MeanStd (Res, ivVal.Ne2, 'iso4Ne');
Res = compute_MeanStd (Res, ivVal.dPs2, 'iso1dPs');
Res = compute_MeanStd (Res, ivVal.dPe2, 'iso2dPe');
Res = compute_MeanStd (Res, ivVal.dNs2, 'iso3dNs');
Res = compute_MeanStd (Res, ivVal.dNe2, 'iso4dNe');

% END OF VVCR_MULTIH
end

%% Auxilliary Functions to simplify computation of final values.

% --- Interpret the return structure coming back from GUIs.
function [Res, Ret] = interpret_str (str, guinam, patfile, ResIn)

Res = ResIn;
Ret = 1;

% GUI was closed using OS
if ~isstruct(str)
    Res = false;
    disp(['calc_vvcr: ' guinam ' closed.']);
    return
end

% The exit button has been pressed
if str.Exit == false
    Res = false;
    disp('calc_vvcr: You chose to exit the analysis');
    disp(['    The file ', patfile, ' was not evaluated.']);

% The discard patient button has been pressed
elseif str.Exit == true
    Res = true;
    disp('calc_vvcr: patient discarded pre-analysis.');

% Everything went well
else
    Ret = 0;
end

end
% --- end interpret_str ---

% --- Compute mean and standard deviation of the named input variable
function [Out] = compute_MeanStd (In, Var, nam)

Out = In;

Out.([nam '_Mean']) = mean(Var);
Out.([nam '_StD'])  = std(Var);

end
% --- end compute_MeanStd ---

% --- Compute mean and standard deviation of VVCR from the named input pressure
function [Out] = compute_VVCR (In, Pes, Pmx, nam)

Out = In;

% Normal Ees/Ea or (Pmx-Pes)/Pes
fieldmean = [nam '_VVCRn_Mean'];
fieldstd  = [nam '_VVCRn_StD'];

Out.(fieldmean) = (Pmx./Pes)-1;
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

% Inverse Ea/Ees or Pes/(Pmx-Pes)
fieldmean = [nam '_VVCRi_Mean'];
fieldstd  = [nam '_VVCRi_StD'];

Out.(fieldmean) = Pes./(Pmx-Pes);
Out.(fieldstd)  = std(Out.(fieldmean));
Out.(fieldmean) = mean(Out.(fieldmean));

end
% --- end compute_VVCR ---

% --- Create blank fields in the event an analysis type was skipped.
function [Res] = create_blank_fields (nam, ResIn, VcorrFlag);

Res = ResIn;

fields = {'_VVCRn', '_VVCRi', '_Pmax'};
for i = 1 : 1 : 3
    Res.([nam fields{i} '_Mean']) = 0;
    Res.([nam fields{i} '_StD'])  = 0;
end

fields = {'_AllDat', '_nFit'};
for i = 1 : 1 : 2
    Res.([nam fields{i}]) = 0;
end

if VcorrFlag
    Res.([nam '_Vcorr']) = 0;
end

end
% --- end create_blank_fields ---
