function [Ret1, Ret2] = data_isoseg (iv1PsIdx, iv1NeIdx, time, timeDoub, ...
                                PresDoub, dPmaxIdx, dPminIdx, GUI, Data, ivIdx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iv%1%PsIdx = ivIdx.Ps1; 
%iv%1%NeIdx = ivIdx.Ne1;
%ti%m%e = Data.Time;
%ti%m%eDoub = Data.Time_D;
%Pr%e%sDoub = Data.Pres_D;
%dP%m%axidx = ivIdx.dPmax;
%dP%m%inidx = ivIdx.dPmin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LEGACY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mysz = length(iv1PsIdx);

% Each row holds the data for a single pressure wave; preallocate all
iv1Pres = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
iv1Time = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
if ~GUI
    iv1PsIdx_Doub = zeros(mysz,1);
    iv1NeIdx_Doub = zeros(mysz,1);
else
    iv1PsIdx_Doub = iv1PsIdx;
    iv1NeIdx_Doub = iv1NeIdx;
end

for i = 1: mysz
    % Positive (1st Isovolumic Section)
    % convert index to index of vector containg 2x data points
    P2 = find(round(timeDoub,3)==round(time(dPmaxIdx(i)),3));
%   [P2, dPmaxIdx(i)]
    if ~GUI
        iv1PsIdx_Doub(i) = find(round(timeDoub,3)==round(time(iv1PsIdx(i)),3));
    end
    % keep in mind these are indices of the time vector, not real time points
    iv1Time(i).PosIso(:,1) = (iv1PsIdx_Doub(i):1:P2)'; 
    for j = 1:length(iv1Time(i).PosIso)
        iv1Pres(i).PosIso(j,1) = PresDoub(iv1Time(i).PosIso(j,1)); % reall pressure points [mmHg]
    end

    %Negative (2nd Isovolumic Section)
    % convert index to index of vector containg 2x data points
    P1 = find(round(timeDoub,3)==round(time(dPminIdx(i)),3));
%   [P1, dPminIdx(i)]
    if ~GUI
        iv1NeIdx_Doub(i) = find(round(timeDoub,3)==round(time(iv1NeIdx(i)),3));
    end
    iv1Time(i).NegIso(:,1) = (P1:1:iv1NeIdx_Doub(i))';
    for j = 1:length(iv1Time(i).NegIso)
        iv1Pres(i).NegIso(j,1) = PresDoub(iv1Time(i).NegIso(j,1));
    end
end

Ret1.P = iv1Pres;
Ret1.T = iv1Time;
if ~GUI
    Ret2.i1d = iv1PsIdx_Doub;
    Ret2.i2d = iv1NeIdx_Doub;
end
