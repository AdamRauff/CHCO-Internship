function [Ret1, Ret2] = data_isoseg (GUI, Data, ivIdx)

% convert time, pressure indices to 2x data points indicies

Ret2 = ivIdx;
mysz = length(ivIdx.Ps1);

% Storing isovolumetric data in structure with two fields: First field is 
% positive isovolmetric, storing all the points that lie on the left side, 
% the positive slope side of the pressure wave, and the second field stores 
% the points on the negative slope side.
%
% Each row holds the data for a single pressure wave; preallocate all
% Note that iv1Time structure holds indices, not actual time points; however
% the iv1Pres structure holds pressures.

Ret1.iv1Time = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
Ret1.iv1Pres = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
Ret1.iv2Time = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));
Ret1.iv2Pres = struct('PosIso', cell(mysz,1), 'NegIso', cell(mysz,1));

% When called from VVCR_MULTI, we build the doubled-up indexes. When called
% from a GUI, we just use those that were already created.
if ~GUI
    Ret2.Ps1_D = zeros(mysz,1);
    Ret2.Ne1_D = zeros(mysz,1);
    Ret2.Pe2_D = zeros(mysz,1);
    Ret2.Ns2_D = zeros(mysz,1);
    Ret2.Ne2_D = zeros(mysz,1);
end

for i = 1: mysz

    % Takaguchi Indices
    % Positive (1st Isovolumic Section)
    P2 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmax(i)),3));
    if ~GUI
        Ret2.Ps1_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ps1(i)),3));
    end

    Ret1.iv1Time(i).PosIso(:,1) = (Ret2.Ps1_D(i):1:P2)'; 
    Ret1.iv1Pres(i).PosIso(:,1) = Data.Pres_D(Ret1.iv1Time(i).PosIso(:,1));

    % Negative (2nd Isovolumic Section)
    P1 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmin(i)),3));
    if ~GUI
        Ret2.Ne1_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ne1(i)),3));
    end

    Ret1.iv1Time(i).NegIso(:,1) = (P1:1:Ret2.Ne1_D(i))';
    Ret1.iv1Pres(i).NegIso(:,1) = Data.Pres_D(Ret1.iv1Time(i).NegIso(:,1));

    % Kind Indicies
    % Positive (1st Isovolumic Section)
    if ~GUI
        Ret2.Pe2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Pe2(i)),3));
    end

    Ret1.iv2Time(i).PosIso(:,1) = (Ret2.Ps1_D(i):1:Ret2.Pe2_D(i))'; 
    Ret1.iv2Pres(i).PosIso(:,1) = Data.Pres_D(Ret1.iv2Time(i).PosIso(:,1));

    % Negative (2nd Isovolumic Section) ***NOTE THIS IS dP/dt***
    if ~GUI
        Ret2.Ns2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Ns2(i)),3));
        Ret2.Ne2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Ne2(i)),3));
    end

    Ret1.iv2Time(i).NegIso(:,1) = (Ret2.Ns2_D(i):1:Ret2.Ne2_D(i))';
    Ret1.iv2Pres(i).NegIso(:,1) = Data.dPdt_D(Ret1.iv2Time(i).NegIso(:,1));

end

