function [ Res, Pat ] = calc_vvcr(PathName, FileName)

% If GUI is true, code will pop up every GUI every time. If GUI is false, it
% will only pop up GUI_GateCheck when the number of maxes & mins is not equal,
% and will not pop up the GUI_Fit* figures. Things go faster, but it's dangerous
% (no checks on the output).
GUI = true;

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
for i = 1:2
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

    % call on GUI.
    if length(Extr.dPmaxVal) == length(Extr.dPminVal)
        if GUI
            GateStr = GUI_GateCheck (PeakStr);
        else
            GateStr = PeakStr;
            GateStr.Exit = 'good';
            GateStr.TotNumWaves = length(Extr.dPmaxVal);
        end
    else
        GateStr = GUI_GateCheck (PeakStr);
    end
        
    [Res, Ret] = interpret_str (GateStr, 'GUI_GateCheck', FileName, Res);
    if Ret
        return;
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
    [ivIdx, ivVal] = data_isoidx (Data_O, Extr);

    % If very few timings were found, filtering may be a problem.
    mingood = min([length(ivIdx.Ps1) length(ivIdx.Ps2)]);
    Found  = double(mingood)/double(Res.A_TotNumWaves);
    if i == 1 && Found < 0.5 && GUI
	    quest = ['Very few cycles gated compared to total number of ' ...
            'cycles. Keep current (filtered) data, load unfiltered data, ' ...
            'or discard patient?'];
        button = questdlg(quest,'Poor Gating','Keep Current', ...
            'Load Unfiltered', 'Discard Patient', 'Load Unfiltered');

        switch button
            case 'Keep Current'
                break;

            case 'Load Unfiltered'
                disp(['calc_vvcr: poor fit gating of max/min cycles. ' ...
                    'Retrying with raw data.']);
                Data_O.FiltPres = Data_O.Pres;
                Data_O.FiltdPdt = Data_O.dPdt;
                Data_O.FiltdP2t = Data_O.dP2t;
                Data_O.Pres = Data_O.OrigPres; 
                Data_O.dPdt = Data_O.OrigdPdt;
                Data_O.dP2t = Data_O.OrigdP2t;
                Data_O = rmfield(Data_O, 'OrigPres');
                Data_O = rmfield(Data_O, 'OrigdPdt');
                Data_O = rmfield(Data_O, 'OrigdP2t');

            case 'Discard Patient'
                Res = true;
                Pat = true;
                disp('calc_vvcr: patient discarded pre-analysis.');
                return
        end

    else
        break;

    end     
end

RunT = logical(length(ivIdx.Ps1));
RunK = logical(length(ivIdx.Ps2));
RunV = logical(length(ivIdx.Ps3));

% if there were no good pressure waveforms left, then skip patient
if ~RunT && ~RunK && ~RunV

    % set all output variables to true, return to runAll
    [Res, Pat] = deal(true);
    return

end

%% (5) Fourier Series interpolation and finding isovolumic segments for fitting
[Data] = data_double (Data_O, ivIdx);

[ivSeg, ivIdx] = data_isoseg (false, Data, ivIdx);

% Get Pes processing done prior to fits: this is for all measured Pes, not just
% the per-cycle (accepted) fits... that might be bad? Not sure. Below, when VVCR
% is computed, only Pes for accepted cycles gets passed.
Res = compute_MeanStd (Res, Data.PesP, 'B_PesP');

% In all the fits below, FHEAD values prior to entry to each RunX section
% define the headers for returned fields. The X_ in front of each name allows
% them to be sorted as desired when returned to runAll. Note that other fields
% already set (like the B_PesP above) also have X_ in front to allow sorting.

%% (6) Perform Takeuchi fit(s), put up check GUI, and compute return quantities
% FitT is "new" (time-normalized) w/new Pes, FitO is same fit (Pmax) with old
% dog Pes.
FHEAD = 'E_TakeNorm';
if RunT
    % Frequency based on individual cycle length works better, that is now done
    % in the fit_ routine on a per-cycle basis.
    ICS.Pres = Data.Pres;
    ICS.dPmaxIdx = ivIdx.dPmax1;
    ICS.dPminIdx = ivIdx.dPmin1;

    % Keep non-normalized fit for now, but don't include its output anymore.
    % When time allows, remove the call below [FitO] and its processing in the
    % GUI.
    [FitT, ivSeg, PlotT] = fit_takeuchi (ivSeg, Data, ICS, 1);

    % Call the Takeuchi Fit Check GUI
    TStr.FitT = FitT;
    if GUI
        TStr.Plot = PlotT;
        RetT = GUI_FitTakeuchi (TStr, Data, ivIdx, ivVal, ivSeg);
        [Res, Ret] = interpret_str (RetT, 'GUI_FitTakeuchi', FileName, Res);
        if Ret
            return;
        end
    else
        RetT = TStr;
    end
   
    BadCycT = RetT.FitT.BadCyc;

    Res.([FHEAD '_AllDat']) = RetT.FitT;
    Res.([FHEAD '_nFit']) = sum(~BadCycT);
    Res.([FHEAD '_Vcorr']) = sum(RetT.FitT.VCyc);

    GOOD_PmxT = RetT.FitT.PIsoMax(BadCycT~=1);
    Res = compute_MeanStd (Res, GOOD_PmxT, [FHEAD '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycT~=1), GOOD_PmxT, FHEAD);
else
    [Res] = create_blank_fields (FHEAD, Res, true);
end

%% (7) Perform Takeuchi fit w/Vanderpool Landmarks, put up check GUI, and 
% compute return quantities. FitV uses "new" ICs (w/fitting limits and new as
% with above. No "old" (unconstrained) fit here.
FHEAD = 'F_Vander';
if RunV
    % Frequency based on individual cycle length works better, that is now done
    % in the fit_ routine on a per-cycle basis.
    ICS.Pres = Data.Pres;
    ICS.dPmaxIdx = ivIdx.dPmax3;
    ICS.dPminIdx = ivIdx.dPmin3;

    [FitV, ivSeg, PlotV] = fit_takeuchi (ivSeg, Data, ICS, 2);

    % Call the Vanderpool Fit Check GUI
    VStr.FitV = FitV;
    if GUI
        VStr.Plot = PlotV;
        RetV = GUI_FitVanderpool (VStr, Data, ivIdx, ivVal, ivSeg);
        [Res, Ret] = interpret_str (RetV, 'GUI_FitVanderpool', FileName, Res);
        if Ret
            return;
        end
    else
        RetV = VStr;
    end

    BadCycV = RetV.FitV.BadCyc;

    Res.([FHEAD '_AllDat']) = RetV.FitV;
    Res.([FHEAD '_nFit']) = sum(~BadCycV);
    Res.([FHEAD '_Vcorr']) = sum(RetV.FitV.VCyc);
    
    GOOD_PmxV = RetV.FitV.PIsoMax(BadCycV~=1);
    Res = compute_MeanStd (Res, GOOD_PmxV, [FHEAD '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycV~=1), GOOD_PmxV, FHEAD);
else
    [Res] = create_blank_fields (FHEAD, Res, true);
end

%% (8) Perform Kind fit, put up check GUI, and compute return quantities
% FitK is unweighted residuals with time normalization. FitN has contraction
% error weighted to be (roughly) the same as relaxtion error.
FHEAD = 'C_KindNorm';
if RunK

    % FitN method can also be eliminated (or converted into experimental method
    % for computing 
    [FitK, PlotK] = fit_kind (ivSeg, ivIdx, Data, 0);

    KStr.FitK = FitK;
    if GUI
        KStr.Plot = PlotK;
        RetK = GUI_FitKind (KStr, Data, ivIdx, ivVal, ivSeg);
        [Res, Ret] = interpret_str (RetK, 'GUI_FitKind', FileName, Res);
        if Ret
            return;
        end
    else
        RetK = KStr;
    end

    BadCycK = RetK.FitK.BadCyc;

    Res.([FHEAD '_AllDat']) = RetK.FitK;
    Res.([FHEAD '_nFit']) = sum(~BadCycK);

    GOOD_PmxK = RetK.FitK.RCoef(BadCycK~=1,1);
    Res = compute_MeanStd (Res, GOOD_PmxK, [FHEAD '_Pmax']);
    Res = compute_VVCR (Res, Data.PesP(BadCycK~=1), GOOD_PmxK, FHEAD);
else
    [Res] = create_blank_fields (FHEAD, Res, false);
end

% Landmark Values
Res = compute_MeanStd (Res, ivVal.dPmax2, 'G_dPmax_P');
Res = compute_MeanStd (Res, ivVal.dPmin2, 'G_dPmin_P');

Res = compute_MeanStd (Res, ivVal.Ps1, 'isoP_T_Ps');
Res = compute_MeanStd (Res, ivVal.Ne1, 'isoP_T_Ne');

Res = compute_MeanStd (Res, ivVal.Ps2, 'isoP_K_Ps');
Res = compute_MeanStd (Res, ivVal.Pe2, 'isoP_K_Pe');
Res = compute_MeanStd (Res, ivVal.Ns2, 'isoP_K_Ns');
Res = compute_MeanStd (Res, ivVal.Ne2, 'isoP_K_Ne');

Res = compute_MeanStd (Res, ivVal.Ps3, 'isoP_V_Ps');
Res = compute_MeanStd (Res, ivVal.Ne3, 'isoP_V_Ne');

% Landmark Timing
norm = Data.time_step./Data.time_end;

Res = compute_MeanStd (Res, norm*ivIdx.dPmax2, 'G_dPmax_T');
Res = compute_MeanStd (Res, norm*ivIdx.dPmin2, 'G_dPmin_T');

Res = compute_MeanStd (Res, norm*ivIdx.Ps1, 'isoT_T_Ps');
Res = compute_MeanStd (Res, norm*ivIdx.Ne1, 'isoT_T_Ne');

Res = compute_MeanStd (Res, norm*ivIdx.Ps2, 'isoT_K_Ps');
Res = compute_MeanStd (Res, norm*ivIdx.Pe2, 'isoT_K_Pe');
Res = compute_MeanStd (Res, norm*ivIdx.Ns2, 'isoT_K_Ns');
Res = compute_MeanStd (Res, norm*ivIdx.Ne2, 'isoT_K_Ne');

Res = compute_MeanStd (Res, norm*ivIdx.Ps3, 'isoT_V_Ps');
Res = compute_MeanStd (Res, norm*ivIdx.Ne3, 'isoT_V_Ne');

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

Res.([nam '_nFit']) = 0;
Res.([nam '_AllDat']).Err = 0;

if VcorrFlag
    Res.([nam '_Vcorr']) = 0;
end

end
% --- end create_blank_fields ---
