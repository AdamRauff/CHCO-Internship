function [Ret1, Ret2] = data_isoseg (GUI, Data, ivIdx)

% convert time, pressure indices to 2x data points indicies
Ret2 = ivIdx;
mysz1 = length(ivIdx.Ps1);
mysz2 = length(ivIdx.Ps2);
mysz3 = length(ivIdx.Ps3);

% Storing isovolumetric data in structure with two fields: First field is 
% positive isovolmetric, storing all the points that lie on the left side, 
% the positive slope side of the pressure wave, and the second field stores 
% the points on the negative slope side.
%
% Each row holds the data for a single pressure wave; preallocate all
% Note that iv1Time structure holds indices, not actual time points; however
% the iv1Pres structure holds pressures.

Ret1.iv1Time = struct('PosIso', cell(mysz1,1), 'NegIso', cell(mysz1,1));
Ret1.iv1Pres = struct('PosIso', cell(mysz1,1), 'NegIso', cell(mysz1,1));
Ret1.iv2Time = struct('PosIso', cell(mysz2,1), 'NegIso', cell(mysz2,1));
Ret1.iv2Pres = struct('PosIso', cell(mysz2,1), 'NegIso', cell(mysz2,1));
Ret1.iv3Time = struct('PosIso', cell(mysz3,1), 'NegIso', cell(mysz3,1));
Ret1.iv3Pres = struct('PosIso', cell(mysz3,1), 'NegIso', cell(mysz3,1));

% When called from VVCR_MULTI, we build the doubled-up indexes. When called
% from a GUI, we just use those that were already created.
if ~GUI
    Ret2.Ps1_D = zeros(mysz1,1);
    Ret2.Ne1_D = zeros(mysz1,1);
    Ret2.Ps2_D = zeros(mysz2,1);
    Ret2.Pe2_D = zeros(mysz2,1);
    Ret2.Ns2_D = zeros(mysz2,1);
    Ret2.Ne2_D = zeros(mysz2,1);
    Ret2.dPmin2_D = zeros(mysz2,1);
    Ret2.Ps3_D = zeros(mysz3,1);
    Ret2.Ne3_D = zeros(mysz3,1);
end

% Takaguchi Indices - every variable below is an index, EXCEPT for the
% iv(x)Pres variables, which are actual values of pressure.
for i = 1: mysz1

    % Positive (1st Isovolumic Section). P2 is an index, not a value, as is
    % Ps1_D, below.
    P2 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmax1(i)),3));
    if ~GUI
        Ret2.Ps1_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ps1(i)),3));
    end

    % iv1Time is indices.
    Ret1.iv1Time(i).PosIso(:,1) = (Ret2.Ps1_D(i):1:P2)'; 
    Ret1.iv1Pres(i).PosIso(:,1) = Data.Pres_D(Ret1.iv1Time(i).PosIso(:,1));

    % Negative (2nd Isovolumic Section). Again here, P1 and Ne1_D(i) are
    % indices.
    P1 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmin1(i)),3));
    if ~GUI
        Ret2.Ne1_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ne1(i)),3));
    end

    % iv1Time is again indices.
    Ret1.iv1Time(i).NegIso(:,1) = (P1:1:Ret2.Ne1_D(i))';
    Ret1.iv1Pres(i).NegIso(:,1) = Data.Pres_D(Ret1.iv1Time(i).NegIso(:,1));
end

% Vanderpool Indices. See above for more explanation. Almost every variable
% below is an index, not a value.
for i = 1: mysz3

    % Positive (1st Isovolumic Section)
    P2 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmax3(i)),3));
    if ~GUI
        Ret2.Ps3_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ps3(i)),3));
    end

    Ret1.iv3Time(i).PosIso(:,1) = (Ret2.Ps3_D(i):1:P2)'; 
    Ret1.iv3Pres(i).PosIso(:,1) = Data.Pres_D(Ret1.iv3Time(i).PosIso(:,1));

    % Negative (2nd Isovolumic Section)
    P1 = find(round(Data.Time_D,3) == round(Data.Time(ivIdx.dPmin3(i)),3));
    if ~GUI
        Ret2.Ne3_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ne3(i)),3));
    end

    Ret1.iv3Time(i).NegIso(:,1) = (P1:1:Ret2.Ne3_D(i))';
    Ret1.iv3Pres(i).NegIso(:,1) = Data.Pres_D(Ret1.iv3Time(i).NegIso(:,1));   
end

% Kind Indicies. See above for more explanation. Almost every variable
% below is an index, not a value.
for i = 1: mysz2

    % Positive (1st Isovolumic Section)
    if ~GUI
        Ret2.Ps2_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.Ps2(i)),3));
        Ret2.Pe2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Pe2(i)),3));
    end

    Ret1.iv2Time(i).PosIso(:,1) = (Ret2.Ps2_D(i):1:Ret2.Pe2_D(i))'; 
    Ret1.iv2Pres(i).PosIso(:,1) = Data.Pres_D(Ret1.iv2Time(i).PosIso(:,1));

    % Negative (2nd Isovolumic Section) ***NOTE THIS IS dP/dt***
    if ~GUI
        Ret2.dPmin2_D(i) = find(round(Data.Time_D,3) == ...
            round(Data.Time(ivIdx.dPmin2(i)),3));
        Ret2.Ns2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Ns2(i)),3));
        Ret2.Ne2_D(i) = find(round(Data.Time_D,3) == ... 
            round(Data.Time(ivIdx.Ne2(i)),3));
    end

    Ret1.iv2Time(i).NegIso(:,1) = (Ret2.Ns2_D(i):1:Ret2.Ne2_D(i))';
    Ret1.iv2Pres(i).NegIso(:,1) = Data.Pres_D(Ret1.iv2Time(i).NegIso(:,1));
    Ret1.iv2dPdt(i).NegIso(:,1) = Data.dPdt_D(Ret1.iv2Time(i).NegIso(:,1));
end

%% NEED CHECKS HERE FOR EMPTY SEGMENTS (JUST IN CASE)...
% Find empty segments and remove them from both the iv(x) vectors as well as
% the iv(XXX) vectors as well. The checks at the end of data_isoidx are now
% functionalized so that they can be called to here remove data from the
% iv(XXX) vectors, along with another subfunction that cleans the iv(X)Pres
% /iv(X)Time vectors as well.