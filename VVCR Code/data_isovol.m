function [RetVal] = isovol_data (Iso1StIdx, Iso2StIdx, time, timeDoub, ...
                                PresDoub, dPmaxIdx, dPminIdx, GUI)

mysz = length(Iso1StIdx);

% Each row holds the data for a single pressure wave; preallocate all
isovolPres = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
isovolTime = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
if ~GUI
    Iso1StIdx_Doub = zeros(mysz,1);
    Iso2StIdx_Doub = zeros(mysz,1);
else
    Iso1StIdx_Doub = Iso1StIdx;
    Iso2StIdx_Doub = Iso2StIdx;
end

for i = 1: mysz
    % Positive (1st Isovolumic Section)
    % convert index to index of vector containg 2x data points
    P2 = find(round(timeDoub,3)==round(time(dPmaxIdx(i)),3));
%   [P2, dPmaxIdx(i)]
    if ~GUI
        Iso1StIdx_Doub(i) = find(round(timeDoub,3)==round(time(Iso1StIdx(i)),3));
    end
    % keep in mind these are indices of the time vector, not real time points
    isovolTime(i).PosIso(:,1) = (Iso1StIdx_Doub(i):1:P2)'; 
    for j = 1:length(isovolTime(i).PosIso)
        isovolPres(i).PosIso(j,1) = PresDoub(isovolTime(i).PosIso(j,1)); % reall pressure points [mmHg]
    end

    %Negative (2nd Isovolumic Section)
    % convert index to index of vector containg 2x data points
    P1 = find(round(timeDoub,3)==round(time(dPminIdx(i)),3));
%   [P1, dPminIdx(i)]
    if ~GUI
        Iso2StIdx_Doub(i) = find(round(timeDoub,3)==round(time(Iso2StIdx(i)),3));
    end
    isovolTime(i).NegIso(:,1) = (P1:1:Iso2StIdx_Doub(i))';
    for j = 1:length(isovolTime(i).NegIso)
        isovolPres(i).NegIso(j,1) = PresDoub(isovolTime(i).NegIso(j,1));
    end
end

RetVal.P = isovolPres;
RetVal.T = isovolTime;
if ~GUI
    RetVal.i1d = Iso1StIdx_Doub;
    RetVal.i2d = Iso2StIdx_Doub;
end
