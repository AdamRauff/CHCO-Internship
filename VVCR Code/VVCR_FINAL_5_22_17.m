function [ AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file, numPeaks, STD_Pes, STD_PMX, TotNumWaves] = VVCR_FINAL_10_10( PathName, FileName)


% determine if file is from calf or humans to apply apprpriate loadp
% function

% if thrid digit/entry of filename is numeric == human file
% otherwise, calf file --> use calf loadp function
if FileName(1) == 'H' && ischar(FileName(2)) && ~isnan(str2double(FileName(3)));
    [Pres, dPdt, ~, Pnam, Pmrn, file, ~, ~]=loadp_10_10_16(PathName,FileName,100);
else
    [Pres, dPdt, Rvals, file, ~]=load_calf_p_5_17_17(PathName,FileName,100);
end
    
% check to see if RV array was NOT found by loaddp
if length(Pres) == 1 && Pres == 0
    
    % set all outputs to true --> skip file
    AVG_Pes = true;
    AVG_Pmax = true;
    VVCR_UT = true;
    VVCR_KH = true;
    Pnam = true;
    Pmrn = true;
    file = true;
    numPeaks = true;
    STD_Pes = true;
    STD_PMX = true;
    TotNumWaves = true;
    
    disp('Loadp did not detect an RV column');
    % return to runAll.m
    return
end
% construct time array (units of 4 milliseconds from catheter machine)
% **************** -------------------------
% This time array works for the human data
% what is the time interval from calf data????
% **************** -------------------------
time_end=0.004*size(Pres,1); time=0.004:0.004:time_end;


%% find dP/dt maxima and minima.

% flip data. used for finding Minima
DataInv = (-1)*dPdt;

% The peaks must exceed 100 [mmHg/s] and be
% separated by the number of total seconds of the sample*2-1.
[pks, pksT] = findpeaks(dPdt, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+6)); % 

% find peaks of inverted data
[~, MinIdx] = findpeaks(DataInv, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+6));

%finding the local minima values in dP/dt
Minima = dPdt(MinIdx); 

% make sure only to analyze complete waveforms:
% data must begin with a dP/dt max, and end with dP/dt min
if MinIdx(1) < pksT(1)
    MinIdx(1) = [];
    Minima(1) = [];
end

if MinIdx(end) < pksT(end)
    pksT(end) = [];
    pks(end) = [];
end

%% GUI for Manual Deletion of Mins and Maxs

% pre-allocate structure
PeakStruct = struct('Data',cell(3,1), 'Min',cell(3,1), 'Max', cell(3,1), 'IM', cell(3,1));

% The GUI allows the user to manually delete inappropriate minima and
% maxima. In order to do that, The following information must be passed
% onto the GUI:

% Pressure and derivative
PeakStruct(1).Data = time;
PeakStruct(2).Data = Pres;
PeakStruct(3).Data = dPdt;

% Minima
PeakStruct(1).Min = MinIdx;
PeakStruct(2).Min = Minima;

% Maxima
PeakStruct(1).Max = pksT;
PeakStruct(2).Max = pks;

% pictures
Green_Check = imread('check.png');
PeakStruct(1).IM = Green_Check;

Red_X = imread('ex.png');
PeakStruct(2).IM = Red_X;

% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
PeakStruct2 = GUI_No_Peaks_10_10(PeakStruct);

clear PeakStruct Green_Check Red_X DataInv

% if the exit button has been pressed
if PeakStruct2(1).Max == false
    
    % set all output variables to false
    AVG_Pes = false;
    AVG_Pmax = false;
    VVCR_UT = false;
    VVCR_KH = false;
    Pnam = false;
    Pmrn = false;
    file = false;
    numPeaks = false;
    STD_Pes = false;
    STD_PMX = false;
    TotNumWaves = false;
    
    disp('You chose to exit the analysis');
    disp(['The file ', FileName, ' was not evaluated!']);
    % return to runAll.m
    return

% if the discard patient button has been pressed
elseif PeakStruct2(3).Max == true
    
     % set all output variables to true
    AVG_Pes = true;
    AVG_Pmax = true;
    VVCR_UT = true;
    VVCR_KH = true;
    Pnam = true;
    Pmrn = true;
    file = true;
    numPeaks = true;
    STD_Pes = true;
    STD_PMX = true;
    TotNumWaves = true;
    
    % return to runAll.m
    return

    % otherwise
else
    
    % update minima and maxima per user filter GUI
    pksT = PeakStruct2(1).Max;
    pks = PeakStruct2(2).Max;

    MinIdx = PeakStruct2(1).Min;
    Minima = PeakStruct2(2).Min;
    
    % obtain number of total waveforms
    TotNumWaves = PeakStruct2(3).Max;
end


%%  Find EDP: 0.2*dP/dt max
    
% Pre allocate variables
EDP = length(pksT);
EDP_T = length(pksT);
EDP_NT = length(pksT); % N for Negative (on the other side of pressure wave, negative derivative)
EDP_N = length(pksT);
bad_curve = [];

% scroll through all maxima
for i = 1:length(pksT)

    EDi = pksT(i);
    
    % while dPdt(EDi) > 0.20*dPdt(pksT(i))
    while dPdt(EDi) > 0.20*pks(i)
        EDi = EDi - 1;
        if EDi == 0 && i == 1
            % the first dP/dt max is too early in the data, and the
            % pressure wave does not contain enough information to
            % include. 

            % Remove the First Maximum and Minimum ****
            EDi = MinIdx(1); % set EDi to be index of a minimum to force exit of while loop
            bad_curve = [bad_curve, i]; % add to list of bad curves
        end
    end

    % sometimes EDi is 1 or 2 time points away from dP/dt max.
    % This is usually to a step-like shape of the pressure curve, because
    % of the error associated with the physcial system of the catheter.
    % if EDi is less then or equal to 3 points away from dP/dt max
    if abs(EDi-pksT(i)) <= 3 
       EDi = EDi - 1; % bump EDi one time point back

       % Continuation mark
       CONT_MARK = true; % this is a flag that keeps track of the 3 point behind EDi
       % try while loop again, additng the condition that it must be more than 4 points away, and the 3 points before EDi must also be below the peak 
       while CONT_MARK
            if dPdt(EDi) <= 0.20*pks(i) && dPdt(EDi - 1) < 0.20*pks(i) && dPdt(EDi - 2) < 0.20*pks(i) && dPdt(EDi - 3) < 0.20*pks(i) && abs(EDi-pksT(i)) > 3
                CONT_MARK = false;
            end
            if dPdt(EDi) <= 0.20*pks(i) && dPdt(EDi - 1) < 0.20*pks(i) && dPdt(EDi - 2) < 0.20*pks(i) && dPdt(EDi - 3) < 0.20*pks(i) && abs(EDi-pksT(i)) <= 3
                bad_curve = [bad_curve, i]; % add to list of bad curves
                break
            end
            EDi = EDi - 1;
       end
    end

    % assign EDP values
    EDP(i) = Pres(EDi);
    EDP_T(i) = EDi;

    % find point on the other side of the pressure wave, with the
    % same height (pressure) as EDP. EDP_N - EDP Negative
    EDP_Ni = EDP_T(i)+15;

    % round values to nearest tenth
    while round(Pres(EDP_Ni),1) > round(EDP(i),1)
        EDP_Ni = EDP_Ni+1;

        if EDP_Ni == length(Pres)
            % the last dPdt min in the data is part of a pressure
            % wave that is not fully contained in the sample

            EDP_Ni = EDP_T(i)-10; % set EDP_Ni to be 10 data points (40 ms) prior to the EDP to force exit while loop
            bad_curve = [bad_curve, i]; % add to list of bad curves
        end
        
        % if algorithm unable to find the negative EDP, and continues to
        % search into proceeding waveforms, cut it off and get rid of that
        % waveform
        
        % first attain sign of dP/dt
        TmpSign = sign(dPdt(EDP_Ni));
        
        % if the derivative is positive, and the time is past dP/dt min (+10 time increments),
        % erase waveform
        if TmpSign == 1 
           if EDP_Ni-MinIdx(i) > 10
               bad_curve = [bad_curve, i]; % add to list of bad curves
               EDP_Ni = EDP_T(i)-5; % set EDP_Ni to be 10 data points (40 ms) prior to the EDP to force exit while loop
           end
        end
    end

    % find which point is closest to the pressure value of EDP
    tempNeg_EDPs = EDP_Ni-3:1:EDP_Ni+3; % create local neighborhood
    
    % if the last pressure wave form is being evaluated
    if i == length(EDP)
        
        % check that the temporary neighborhood does not exceed the size of
        % the pressure / time vector
        while max(tempNeg_EDPs)>length(Pres)
            tempNeg_EDPs(end) = [];
        end
    end
    
    % calculate the difference between the local neighborhood
    % (tempNeg_EDPs) and EDP
    EDP_Diffs = abs(EDP(i)-Pres(tempNeg_EDPs));

    % find minimum difference
    [~, tempInds] = min(EDP_Diffs);

    % assign Negative EDP values
    EDP_NT(i) = tempNeg_EDPs(tempInds); % keep in mind this is an indicie of the time vector
    EDP_N(i) = Pres(tempNeg_EDPs(tempInds)); % pressure in mmHg

    % --------------------------------------------------------------------------
    % if the EDP_NT < MinIdx(i) or EDP_NT is only 4 point ahead of the min, then we got a problem. no
    % isovolumic points on the negative side of the curve. 
    if EDP_NT(i) <= MinIdx(i) || (EDP_NT(i) > MinIdx(i) && abs(EDP_NT(i)-MinIdx(i)) <= 3)
        disp(['Curve # ',num2str(i), ' has smaller EDP_NT then the min!!!']);
        % get rid of curve if it is not already marked
        if isempty(find(bad_curve==i,1))
            bad_curve = [bad_curve, i];
        end
        % ask Hunter about this scenario
        % see patient HA000251.&05 (2nd waveform) for example
    end
    % ----------------------------------------------------------------------------

end

tempNum = length(bad_curve);

% make sure bad_curve has no redundencies
bad_curve = sort(unique(bad_curve));

% Remove bad curves
if ~isempty(bad_curve)
    % loop through bad_curves vector
    for i = tempNum:-1:1
        EDP(bad_curve(i)) = [];
        EDP_T(bad_curve(i)) = [];
        EDP_NT(bad_curve(i)) = [];
        EDP_N(bad_curve(i)) = [];
        pksT(bad_curve(i)) = [];
        pks(bad_curve(i)) = [];
        MinIdx(bad_curve(i)) = [];
        Minima(bad_curve(i)) = [];
    end
end

% if there were no good pressure waveforms left, then skip patient
if isempty(EDP)
    
     % set all output variables to true
    AVG_Pes = true;
    AVG_Pmax = true;
    VVCR_UT = true;
    VVCR_KH = true;
    Pnam = true;
    Pmrn = true;
    file = true;
    numPeaks = true;
    STD_Pes = true;
    STD_PMX = true;
    TotNumWaves = true;
    
    % return to runAll.m
    return

    % otherwise
end
%% FOURIER SERIES INTERPOLATION
% 
% increase the number of data points to provide a tighter fit of the 
% sinusoid that is fitted to the isovolumic points
PresDoub = interpft(Pres,length(Pres)*2);
timeDoub = 0.002:0.002:time(end);

%% Pes = 30 ms prior to dp/dt min %%%%%%%%%%%%%
dtmin_30=time(MinIdx)-0.03; %30 ms prior to dp/dt min

% Using the interpolated data - an exact point can be reached
P_esTimes = timeDoub(uint16(dtmin_30/0.002)); % convert to integer (can't use double as index). Obtained from interpolated data
P_es = PresDoub(uint16(dtmin_30/0.002));

%% Obtaining Isovolumetric points

% storing isovolumetric data in structure with two fields:
% First field is positive isovolumetric, storing all the points that lie on
% the left side, the positive slope side of the pressure wave, and the
% second field stores the points on the negative slope side.

% Each row holds the data for a single pressure wave
isovol = struct('PosIso',cell(length(EDP),1),'NegIso',cell(length(EDP),1));
isovoltime = struct('PosIso',cell(length(EDP),1),'NegIso',cell(length(EDP),1));

% pre-allocate
EDP_Doub = zeros(length(EDP),1);
EDP_NT_Doub = zeros(length(EDP),1);

for i = 1: length(EDP)
    % Positive
    % convert index to index of vector containg 2x data points
    EDP_Doub(i) = find(round(timeDoub,3)==round(time(EDP_T(i)),3)); 
    P2 = find(round(timeDoub,3)==round(time(pksT(i)),3));
    isovoltime(i).PosIso(:,1) = (EDP_Doub(i):1:P2)'; % keep in mind these are indices of the time vector, not real time points
    for j = 1:length(isovoltime(i).PosIso)
        isovol(i).PosIso(j,1) = PresDoub(isovoltime(i).PosIso(j,1)); % reall pressure points [mmHg]
    end
    
    %Negative
    % convert index to index of vector containg 2x data points
    P1 = find(round(timeDoub,3)==round(time(MinIdx(i)),3));
    EDP_NT_Doub(i) = find(round(timeDoub,3)==round(time(EDP_NT(i)),3));
    isovoltime(i).NegIso(:,1) = (P1:1:EDP_NT_Doub(i))';
    for j = 1:length(isovoltime(i).NegIso)
        isovol(i).NegIso(j,1) = PresDoub(isovoltime(i).NegIso(j,1));
    end
end

%% SINUSOIDAL FITTING
%%%%%%%finally, fitting sinusoid to isovolumetric regions%%%%%%%%%

% pre - allocate 
c_tot2 = zeros(length(EDP),4);
P_max2 = zeros(length(EDP),1);
waveFit = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);

totIsoTimePoints = [];
totIsoPresPoints = [];

% freq is an initial condition that does not rquire individual wave
% calculation --> executed outside loop
% frequnecy is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
Freq = double(((2*pi)*TotNumWaves)/(time_end));

% scroll through the number of rows (pressure waves) in the
% structures: isovoltime and isovol
for i = 1:length(EDP)
    
    WaveTs = [timeDoub(isovoltime(i).PosIso)'; timeDoub(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];
    
    % this equation is from Naeiji et al, single beat method of VVC
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 
    
    % deriving the initial values from the data
    % mean - average pressure value between dp/dt max and min (top of
    % curve)
    T1 = pksT(i);
    T2 = MinIdx(i);
    Mea = mean(double(Pres(T1:T2)));
    
    % Amplitude is twice the mean
    Amp = double(1.8*Mea);
    
    % keep in mind this means the initial conditions of every wave fit may
    % be slightly different, While values entered via GUI make ICs same for
    % all waves.
    c2=[Mea, Amp, Freq, -0.5];

    [c,resnorm,~]=lsqnonlin(sin_fun2,c2); %least squares fitting
    
    Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));
    
    % r^2 value
    r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;
    
    % if the fit of the wave was bad, mark that wave
    if r_square2(i) <0.90
       waveFit(i) = 1; 
    end
    
    %getting all the c values in a matrix
    c_tot2(i,:)=c; 
    
    %first equation pmax, A+B
    P_max2(i)=c(1)+abs(c(2)); 
   
    % AR 6/5/17 -----------------------------------------------
    % adding points succesively to beginning of systole to make better fit
    % of sick patients with wide curves
    
    % obtain maximum pressure point on actual curve
    PresMax = max(PresDoub(isovoltime(i).PosIso(1,1):1:isovoltime(i).NegIso(end,1)));
    if r_square2(i) > 0.80 && P_max2(i) < PresMax
       
        % keep count of how many points added to systole side
        count = 0;
        while P_max2(i) < PresMax
        
            % add point to isovoltime(i).PosIso and corresponding isovol(i).PosIso
            isovoltime(i).PosIso = [(isovoltime(i).PosIso(1,1))-1, isovoltime(i).PosIso];
            isovol(i).PosIso = [PresDoub(isovoltime(i).PosIso(1,1)), isovol(i).PosIso];

            % update Wave(x)s variables
            WaveTs = [timeDoub(isovoltime(i).PosIso)'; timeDoub(isovoltime(i).NegIso)'];
            WavePs = [isovol(i).PosIso; isovol(i).NegIso];

            % re-fit sinusiod
            % equation from Naeiji et al, single beat method of VVC
            sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

            %least squares fitting
            [c,resnorm,~]=lsqnonlin(sin_fun2,c2); 

            Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

            % r^2 value
            r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

            % if the fit of the wave was bad, mark that wave
            if r_square2(i) <0.90
               waveFit(i) = 1;
            else
                waveFit(i) = 0;
            end
            
            %getting all the c values in a matrix
            c_tot2(i,:)=c; 
    
            %first equation pmax, A+B
            P_max2(i)=c(1)+abs(c(2));
            
            % increment count to keep track of added points
            count = count +1;
            
            % Do not let program add more than 10 points
            if count >= 10 && (P_max2(i) < PresMax || waveFit(i) == 1)
                waveFit(i) = 1;
                disp('Added nine points on systolic side of curve, and Pmax remains short of actual pressure');
                disp(['Wave: ',num2str(i), 'is excluded']);
                break
            end
        end
    end
    
    % ---------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 (last pressure waveform) for an example 
    
    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which can make a good fit (r^2 > 0.99).
    % --------------------------------------------------------------------

    % store the time points and pressure points in one array for easy
    % plotting 
    totIsoTimePoints = [totIsoTimePoints; WaveTs];
    totIsoPresPoints = [totIsoPresPoints; WavePs];
end
%% GUI to change ICs

% pre-allocate structure
PeakStruct = struct('Data',cell(2,1), 'ivt',cell(2,1), 'iv', cell(2,1), 'isoPts', cell(2,1), 'Cs', cell(2,1), 'Misc', cell(2,1), 'EDPs', cell(2,1), 'Crit', cell(2,1));

% Pressure and derivative
PeakStruct(1).Data = timeDoub;
PeakStruct(2).Data = PresDoub;

% isovolumic time
PeakStruct(1).ivt = isovoltime; % times of all isovolumic points\
PeakStruct(2).ivt = time; % passing the old time vector with 1/2 the points. This is used for the buttondownFcn

% isovolumic pressures
PeakStruct(1).iv = isovol;
PeakStruct(2).iv = waveFit; % this keeps track of which waveforms had a bad fit

% passing the time points in one array for ease of plotting
PeakStruct(1).isoPts = totIsoTimePoints;
PeakStruct(2).isoPts = totIsoPresPoints;

PeakStruct(1).Cs = c2; % intial conditions that were first used
PeakStruct(2).Cs = c_tot2; % regression constants from first fit

PeakStruct(1).Misc = EDP; % used as reference for number of peaks
PeakStruct(2).Misc = P_max2; % Pmax values obtained from fit

PeakStruct(1).EDPs = EDP_Doub; % give the time EDP occured - used for buttondownFcn
PeakStruct(2).EDPs = EDP_NT_Doub; % give the time negative EDP occured - used for buttondownFcn

PeakStruct(1).Crit = pksT; % pass the times of the peaks. used in buttownDownFcn. Recall these are indexs of the old time vector
PeakStruct(2).Crit = MinIdx;

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(waveFit==1); % find indices of the bad waves
if ~isempty(indX)
    disp('The following waves did NOT have a good fit (will not be included)');
    disp(['Wave(s): ', num2str(indX')]);
else
    disp('All waves seemed to fit well!');
end

% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
PeakStruct2 = GUI_SINU_FIT_10_03(PeakStruct);

clear PeakStruct

% if the exit button has been pressed
if size(PeakStruct2(1).output,1) == 1 && PeakStruct2(2).output == false
    
    % set all output variables to false
    AVG_Pes = false;
    AVG_Pmax = false;
    VVCR_UT = false;
    VVCR_KH = false;
    Pnam = false;
    Pmrn = false;
    file = false;
    numPeaks = false;
    STD_Pes = false;
    STD_PMX = false;
    
    % return to runAll.m
    return

% if the discard patient button has been pressed
elseif size(PeakStruct2(1).output,1) == 1 && PeakStruct2(2).output == true
    
     % set all output variables to true
    AVG_Pes = true;
    AVG_Pmax = true;
    VVCR_UT = true;
    VVCR_KH = true;
    Pnam = true;
    Pmrn = true;
    file = true;
    numPeaks = true;
    STD_Pes = true;
    STD_PMX = true;
    
    % return to runAll.m
    return

% otherwise
else
    % extract Pmax and the list of well fitted curves from the structure
    % that is returned from GUI
    P_max2 = PeakStruct2(2).output;
    waveFit = PeakStruct2(1).output;
    
    % calculate the number of peaks that were evaluated
    numPeaks = 0;
    for i = 1:length(waveFit)
        if waveFit(i) ~= 1
            numPeaks = numPeaks + 1;
        end
    end
    %OKAY! here are the final values and we can FINALLY calculate VVCR.

    % average Pes for the waves that fit well
    AVG_Pes=mean(P_es(waveFit~=1)); 
    
    % standard deviation of Pes
    STD_Pes = std(P_es(waveFit~=1));

    % average P_max for the waves that fit well
    AVG_Pmax=mean(P_max2(waveFit~=1)); 
    
    % standard deviation of Pmax
    STD_PMX = std(P_max2(waveFit~=1));

    %from Uyen Truongs VVCR paper
    VVCR_UT=AVG_Pes/(AVG_Pmax-AVG_Pes); 

    % Hunter's 
    VVCR_KH=(AVG_Pmax/AVG_Pes)-1;
end

end

