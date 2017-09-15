%% VVCR code
% Adam Rauff
% September 9th, 2016

% Prepare workspace
clc; close all; clear all

% FLAGS: boolean variables used as flags for ease of debugging
FOUR_DOUB_SAMP = true; % plot the fourier interpolated data that samples the pressure signal at 2x frequency = twice the data points
ORIG_PRES_DPDT = true; % plot pressure vs. time and dP/dt vs. time with Minima, Maxima, and +/- EDP
ES_PTS = false; % Plot End systolic points on top of pressure curve
GRAPH_ISOVOL_POINTS = false; % plot the isovolumetric points of each pressure wave
ISOVOL_REG = true; % display the filtered isovolumetric region, from which the sinusoid is fitted for

% obtain text file that contains the data
[FileName,PathName] = uigetfile('*.*');

% apply loadp function to load the data into matlab variables
[Pres, dPdt, Rvals, name, mrn, marks]=loadp(PathName,FileName,100);

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

% clear variables no longer used
clear PeakStruct Green_Check Red_X DataInv

% update minima and maxima per user filter GUI
pksT = PeakStruct2(1).Max;
pks = PeakStruct2(2).Max;

MinIdx = PeakStruct2(1).Min;
Minima = PeakStruct2(2).Min;

%% FOURIER SERIES INTERPOLATION
% 
% increase the number of data points to provide a tighter fit of the 
% sinusoid that is fitted to the isovolumic points
PresDoub = interpft(Pres,length(Pres)*2);
timeDoub = 0.002:0.002:time(end);

% calculate dPdt from interpolated points
dPdt2 = diff(PresDoub)/0.002;

% cut out first and last 0.2 seconds. too noisy from ringing
dPdt2(1:99) = []; 
dPdt2(end-99:end) = []; 

PresDoub_Short = PresDoub(100:end-101);

% create a matching time vector
dPTime = timeDoub(100:end-101);

% asses goodness of fit of fourier interpolation
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
    
    figure, hold on;
    plot(time, dPdt, 'b', dPTime,dPdt2, 'k');
    set(gca,'fontsize',14);
    title('dP/dt Fourier Interpolation','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    legend('Original', 'Fourier Sampled');
    box on
    grid on
    hold off;
end

%%  Find EDP: 0.2*dP/dt max
    
% Pre allocate variables
EDP = length(pksT);
EDP_T = length(pksT);
EDP_NT = length(pksT); % N for Negative (on the other side of pressure wave, negative derivative)
EDP_N = length(pksT);
BAD_First_WAVE = false;
BAD_Last_WAVE = false;

% scroll through all maxima
for i = 1:length(pksT)
    
    % if first or last waveform, use dPdt, otherwise use dPdt2. reasoning
    % is because of noise
    if i == 1 || i == length(pksT) 
        EDi = pksT(i);
        while dPdt(EDi) > 0.20*pks(i)
            EDi = EDi - 1;
            if EDi == 0 && i == 1
                % the first dP/dt max is too early in the data, and the
                % pressure wave does not contain enough information to
                % include. 

                % Remove the First Maximum and Minimum ****
                EDi = MinIdx(1); % set EDi to be index of a minimum to force exit of while loop
                BAD_First_WAVE = true; % Record flag that denotes a bad first wave
            end
        end
    
    % in the case of a bad first wave, just put zeros for now. removed
    % after for loop
    if i == 1 && BAD_First_WAVE == true
        EDP(1) = 0;
        EDP_T(1) = 0;
        EDP_NT(1) = 0;
        EDP_N(1) = 0;
    else
        % assign EDP values
        EDP(i) = Pres(EDi);
        EDP_T(i) = EDi; % This is an index of the time vector

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
                BAD_Last_WAVE = true; % record flag that remarks a bad last wave
            end
        end

        if BAD_Last_WAVE == false
            % find which point is closest to the pressure value of EDP
            tempNeg_EDPs = EDP_Ni-3:1:EDP_Ni+3; % create local neighborhood

            % calculate the difference between the local neighborhood
            % (tempNeg_EDPs) and EDP
            EDP_Diffs = abs(EDP(i)-Pres(tempNeg_EDPs));

            % find minimum difference
            [~, tempInds] = min(EDP_Diffs);

            % assign Negative EDP values
            EDP_NT(i) = tempNeg_EDPs(tempInds); % keep in mind this is an indicie of the time vector
            EDP_N(i) = Pres(tempNeg_EDPs(tempInds)); % pressure in mmHg
            
            % convert the indices into indices of timeDoub
            indX = find(round(timeDoub,3)==round(time(EDP_NT(i)),3)); 
            EDP_NT(i) = indX;

            % convert the indices into indices of timeDoub
            indX = find(round(timeDoub,3)==round(time(EDP_T(i)),3)); 
            EDP_T(i) = indX;
        end
        
        
    end
    
    else % use the pressure and dPdt with double the points
        
        % find time index where the max occurs in dPdt2
        EDi = find(round(dPTime,3)==round(time(pksT(i))),3);
        
        % sometimes the find command doesn't work :(
        if isempty(EDi)
           EDi = uint16(time(pksT(i))/0.002)-105;
           while round(dPTime(EDi),3) ~= round(time(pksT(i)),3)
               EDi = EDi + 1;
           end
        end
        
        while dPdt2(EDi) > 0.20*pks(i)
            EDi = EDi - 1;
        end
        
        % assign EDP values
        EDP(i) = PresDoub_Short(EDi);
        EDP_T(i) = EDi; % this is an index of the dPTime vector

        % find point on the other side of the pressure wave, with the
        % same height (pressure) as EDP. EDP_N - EDP Negative
        EDP_Ni = EDP_T(i)+30;

        % round values to nearest tenth
        while round(PresDoub_Short(EDP_Ni),1) > round(EDP(i),1)
            EDP_Ni = EDP_Ni+1;
        end
        
        % find which point is closest to the pressure value of EDP
        tempNeg_EDPs = EDP_Ni-3:1:EDP_Ni+3; % create local neighborhood

        % calculate the difference between the local neighborhood
        % (tempNeg_EDPs) and EDP
        EDP_Diffs = abs(EDP(i)-PresDoub_Short(tempNeg_EDPs));

        % find minimum difference
        [~, tempInds] = min(EDP_Diffs);

        % assign Negative EDP values
        EDP_NT(i) = tempNeg_EDPs(tempInds); % keep in mind this is an indicie of the dPTime vector
        EDP_N(i) = PresDoub_Short(tempNeg_EDPs(tempInds)); % pressure in mmHg
        
        % convert the indices into indices of timeDoub
        indX = find(round(timeDoub,3)==round(dPTime(EDP_NT(i)),3)); 
        EDP_NT(i) = indX;

        % convert the indices into indices of timeDoub
        indX = find(round(timeDoub,3)==round(dPTime(EDP_T(i)),3)); 
        EDP_T(i) = indX;
    end
    
end

% Remove first Minima and Maxima, as well as first EDP and EDP_T
if BAD_First_WAVE == true
    EDP(1) = [];
    EDP_T(1) = [];
    EDP_NT(1) = [];
    EDP_N(1) = [];
    pksT(1) = [];
    pks(1) = [];
    MinIdx(1) = [];
    Minima(1) = [];
end

% Remove last Minima and Maxima, as well as last EDP and EDP_T
if BAD_Last_WAVE == true
    EDP(end) = [];
    EDP_T(end) = [];
    pksT(end) = [];
    pks(end) = [];
    MinIdx(end) = [];
    Minima(end) = [];
end

% Plot data with Minima, Maxima, EDP, and Negative EDP
if ORIG_PRES_DPDT == true
    
    figure, hold on;
    plot(timeDoub,PresDoub,'b',time(pksT), Pres(pksT), 'ro', time(MinIdx), Pres(MinIdx), 'ko', timeDoub(EDP_T), EDP, 'mo', timeDoub(EDP_NT), EDP_N, 'co'); hold on;
    set(gca,'fontsize',14);
    title('Original Pressure Vs. Time','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'EDP', 'Negative EDP');
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
    legend('Pressure', 'Pes points', 'Interpolated points', 'dP/dt Min');
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

for i = 1: length(EDP)
    % Positive
    isovoltime(i).PosIso(:,1) = (EDP_T(i):1:pksT(i))'; % keep in mind these are indices of the time vector, not real time points
    for j = 1:length(isovoltime(i).PosIso)
        isovol(i).PosIso(j,1) = Pres(isovoltime(i).PosIso(j,1)); % reall pressure points [mmHg]
    end
    
    %Negative
    isovoltime(i).NegIso(:,1) = (MinIdx(i):1:EDP_NT(i))';
    for j = 1:length(isovoltime(i).NegIso)
        isovol(i).NegIso(j,1) = Pres(isovoltime(i).NegIso(j,1));
    end
end

if GRAPH_ISOVOL_POINTS == true
    figure, hold on;
    plot(time, Pres, 'k'); hold on;
    for i = 1:length(EDP)
        plot(time(isovoltime(i).PosIso), isovol(i).PosIso, 'bo',time(isovoltime(i).NegIso), isovol(i).NegIso, 'ro');
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
P_max2 = zeros(length(EDP),1);
waveFit = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);

totIsoTimePoints = [];
totIsoPresPoints = [];

% scroll through the number of rows (pressure waves) in the
% structures: isovoltime and isovol
for i = 1:length(EDP)
    
    WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];
    
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; %this
    % equation is from Naeiji et al, single beat method of VVC
    c2=[15, 215, 10, -1];

    [c,resnorm,residual]=lsqnonlin(sin_fun2,c2); %least squares fitting
    
    Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));
    
    % r^2 value
    r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;
    
    % if the fit of the wave was bad, mark that wave
    if r_square2(i) <0.95
       waveFit(i) = 1; 
    end
    
    c_tot2(i,:)=c; %getting all the c values in a matrix
    
    P_max2(i)=c(1)+2*abs(c(2)); %first equation pmax, A+2B
    
      % ------------------------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 (last pressure waveform) for an example of why
    
    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which makes a good fit (r^2 > 0.99).
    % ---------------------------------------------------------------------------------
    

    % store the time points and pressure points in one array for easy
    % plotting 
    totIsoTimePoints = [totIsoTimePoints; WaveTs];
    totIsoPresPoints = [totIsoPresPoints; WavePs];
end

%% GUI to change ICs

% pre-allocate structure
PeakStruct = struct('Data',cell(2,1), 'ivt',cell(2,1), 'iv', cell(2,1), 'isoPts', cell(2,1), 'Cs', cell(2,1), 'Misc', cell(2,1));

% Pressure and derivative
PeakStruct(1).Data = time;
PeakStruct(2).Data = Pres;

% isovolumic time
PeakStruct(1).ivt = isovoltime; % times of all isovolumic points\

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

% call on GUI. Notice the structure we just made is passed to the GUI, and
% the GUI passes back a refined structure
PeakStruct2 = GUI_SINU_FIT(PeakStruct);

clear PeakStruct
P_max2 = PeakStruct2(2).output;
waveFit = PeakStruct2(1).output;

%OKAY! here are the final values and we can FINALLY calculate VVCR.

% average Pes
AVG_Pes=mean(P_es); 

% average P_max for the waves that fit well
AVG_Pmax=mean(P_max2(waveFit~=1)); 

%this is the eqation from Uyen Truongs VVCR paper
VVCR_UT=AVG_Pes/(AVG_Pmax-AVG_Pes); 

%this is kendalls
VVCRinvert=(AVG_Pmax/AVG_Pes)-1;
