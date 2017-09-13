%% VVCR code
% Adam Rauff
% September 17th, 2016

% Prepare workspace
clc; close all; % clear all

% FLAGS: boolean variables used as flags for ease of debugging
ORIG_PRES_DPDT = false; % plot pressure vs. time and dP/dt vs. time with Minima, Maxima, and +/- EDP
FOUR_DOUB_SAMP = false; % plot the fourier interpolated data that samples the pressure signal at 2x frequency = twice the data points
ES_PTS = false; % Plot End systolic points on top of pressure curve
GRAPH_ISOVOL_POINTS = false; % plot the isovolumetric points of each pressure wave
SINU_FIT = false; % display sinusoidal fit in its own figure after Sinu fit GUI

% obtain text file that contains the data
[FileName,PathName] = uigetfile('*.*');

% apply loadp function to load the data into matlab variables
[Pres, dPdt, rvals, Pnam, Pmrn, file, npath, marks] = loadp(PathName,FileName,100);

% construct time array (units of 4 milliseconds from catheter machine)
time_end=0.004*size(Pres,1); time=0.004:0.004:time_end;

%% find dP/dt maxima and minima.

% flip data. used for finding Minima
DataInv = (-1)*dPdt;

% The peaks must exceed 100 [mmHg/s] and be
% separated by the number of total seconds of the sample*2-1.
[pks, pksT] = findpeaks(dPdt, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+4)); % 

% find peaks of inverted data
[~, MinIdx] = findpeaks(DataInv, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2+4));

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
PeakStruct2 = GUI_No_Peaks(PeakStruct);

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
    
    disp('You chose to exit the analysis');
    disp(['The file ', FileName, ' was not evaluated!']);
    % return to runAll.m
%     return

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
    
    % return to runAll.m
%     return

    % otherwise
else
    % update minima and maxima per user filter GUI
    pksT = PeakStruct2(1).Max;
    pks = PeakStruct2(2).Max;

    MinIdx = PeakStruct2(1).Min;
    Minima = PeakStruct2(2).Min;
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
    % if EDi is less then or equal to 4 points away from dP/dt max
    if abs(EDi-pksT(i)) <= 4 
       EDi = EDi - 1; % bump EDi one time point back

       % Continuation mark
       CONT_MARK = true; % this is a flag that keeps track of the 3 point behind EDi
       % try while loop again, additng the condition that it must be more than 4 points away, and the 3 points before EDi must also be below the peak 
       while CONT_MARK
            if dPdt(EDi) <= 0.20*pks(i) && dPdt(EDi - 1) < 0.20*pks(i) && dPdt(EDi - 2) < 0.20*pks(i) && dPdt(EDi - 3) < 0.20*pks(i) && abs(EDi-pksT(i)) > 4
                CONT_MARK = false;
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
        
        % if the derivative is positive, and the time is past dP/dt min (+4 increments),
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
    if EDP_NT(i) <= MinIdx(i) || (EDP_NT(i) > MinIdx(i) && abs(EDP_NT(i)-MinIdx(i)) <= 4)
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
    
    % return to runAll.m
%     return

    % otherwise
end
% Plot data with Minima, Maxima, EDP, and Negative EDP
if ORIG_PRES_DPDT == true
    
    figure, hold on;
    plot(time,Pres,'b',time(pksT), Pres(pksT), 'ro', time(MinIdx), Pres(MinIdx), 'ko', time(EDP_T), EDP, 'mo', time(EDP_NT), EDP_N, 'co'); hold on;
    set(gca,'fontsize',14);
    title('Original Pressure Vs. Time','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'EDP', 'Negative EDP');
    box on
    grid on
    hold off;

    % Plot derivative of original data
    figure, hold on;
    plot(time,dPdt, 'b', time(pksT), pks, 'ro', time(MinIdx), Minima, 'ko', time(EDP_T), dPdt(EDP_T), 'mo', time(EDP_NT), dPdt(EDP_NT), 'co'); hold on;
    set(gca,'fontsize',14);
    title(' dP/dt Vs. Time','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    legend('dP/dt','Maxima', 'Minima', 'EDP', 'Negative EDP');
    box on
    grid on
    hold off;
end

%% FOURIER SERIES INTERPOLATION
% 
% increase the number of data points to provide a tighter fit of the 
% sinusoid that is fitted to the isovolumic points
PresDoub = interpft(Pres,length(Pres)*2);
timeDoub = 0.002:0.002:time(end);

% assess goodness of fit of fourier interpolation
if FOUR_DOUB_SAMP == true
    figure, hold on;
    plot(time,Pres, 'b', timeDoub, PresDoub, 'k');
    set(gca,'fontsize',14);
    title('Fourier Interpolation','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressure [mmHg]','FontSize',18);
    legend('Original', 'Fourier');
    box on
    grid on
    hold off;
end

%% Pes = 30 ms prior to dp/dt min %%%%%%%%%%%%%
dtmin_30=time(MinIdx)-0.03; %30 ms prior to dp/dt min

% Using the interpolated data - an exact point can be reached
P_esTimes = timeDoub(uint16(dtmin_30/0.002)); % convert to integer (can't use double as index). Obtained from interpolated data
P_es = PresDoub(uint16(dtmin_30/0.002));

if ES_PTS == true
    figure, hold on;
    plot(time,Pres,'b',P_esTimes,P_es,'ro', time(MinIdx), Pres(MinIdx), 'ko');
    set(gca,'fontsize',14);
    title('Pes 30 [ms] prior to Min dP/dt', 'FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressure [mmHg]','FontSize',18);
    legend('Pressure', 'Pes points', 'dP/dt Min');
    box on;
    grid on;
    hold off;
end


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

if GRAPH_ISOVOL_POINTS == true
    figure, hold on;
    plot(timeDoub, PresDoub, 'k'); hold on;
    for i = 1:length(EDP)
        plot(timeDoub(isovoltime(i).PosIso), isovol(i).PosIso, 'bo',timeDoub(isovoltime(i).NegIso), isovol(i).NegIso, 'ro');
    end
    set(gca,'fontsize',14);
    title('Isovolumetric Region','FontSize',20);
    ylabel('Pressure [mmHg]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    legend('Pressure','Positive Isovolumic','Negative IsoVolumic');
    grid on;
    box on;
end

%% SINUSOIDAL FITTING
%%%%%%%finally, fitting sinusoid to isovolumetric regions%%%%%%%%%

% pre - allocate 
c_tot2 = zeros(length(EDP),4);
waveFit = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);
P_max2 = zeros(length(EDP),1);

totIsoTimePoints = [];
totIsoPresPoints = [];

% scroll through the number of rows (pressure waves) in the
% structures: isovoltime and isovol
for i = 1:length(EDP)
    
    WaveTs = [timeDoub(isovoltime(i).PosIso)'; timeDoub(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];
    
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; %this
    % equation is from Naeiji et al, single beat method of VVC
    c2=[15, 215, 10, -1];

    [c,resnorm,residual]=lsqnonlin(sin_fun2,c2); %least squares fitting
    
    Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));
    
    % r^2 value
    r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;
    
    % if the fit of the wave was bad, mark that wave
    if r_square2(i) <0.85
       waveFit(i) = 1; 
    end
    
    c_tot2(i,:)=c; %getting all the c values in a matrix
    
    P_max2(i)=c(1)+abs(c(2)); %first equation pmax, A+B
    
      % ------------------------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 (last pressure waveform) for an example 
    
    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which can make a good fit (r^2 > 0.99).
    % ---------------------------------------------------------------------------------

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
if size(PeakStruct2(1).output,1) == 1 && PeakStruct2(1).output == false
    
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
elseif size(PeakStruct2(1).output,1) == 1 && PeakStruct2(1).output == true
    
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
    STD_PMX = std(P_es(waveFit~=1));

    %from Uyen Truongs VVCR paper
    VVCR_UT=AVG_Pes/(AVG_Pmax-AVG_Pes); 

    % Hunter's 
    VVCR_KH=(AVG_Pmax/AVG_Pes)-1;

    % make plot of sinusoidal fits
    if SINU_FIT == true
        figure, hold on;
        plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints, 'ro');
        set(gca,'fontsize',14);
        title('Sinusoidal Fitting', 'FontSize',20);
        xlabel('Time [s]','FontSize',18);
        ylabel('Pressure [mmHg]','FontSize',18);
        hold on

        % Attain the sinusoid fit for all points (so Pmax can be visualized
        for i = 1:length(EDP)

            % obtain the range of time of each peak
            interval = timeDoub(isovoltime(i).PosIso(1,1)):0.002:timeDoub(isovoltime(i).NegIso(end,1));

            % plug into Naeiji equation that was just solved for 
            FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));

            % find time point corresponding to Pmax
            [~, Idx] = min(abs(FitSinePres-P_max2(i)));

            PmaxT(i) = interval(Idx);

            plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
            hold on;
        end


        % check the range of pressure values of Pmax. if the max p_max value is
        % over 450, rescale y axis to (0, 300), so individual waveforms can be seen
        if max(P_max2) > 450
            ylim([0, 300]);
        end

        legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit', 'Pmax', 'Location','southoutside', 'Orientation', 'horizontal');
        box on;
        grid on;
        hold off;
    end
end