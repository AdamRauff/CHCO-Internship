%% Kendalls old kid  
clc; close all; clear all

% FLAGS: boolean variables used as flags for ease of debugging
ORIG_PRES_DPDT = true; % original pressure vs. time plot and original dP/dt vs. time plot
FOUR_PRES = false; % fourier interpolation of pressure curve
ES_PTS = false; % Plot End systolic points on top of pressure curve
GRAPH_ISOVOL_POINTS = true; % plot the isovolumetric points of each pressure wave
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
[pks, pksT] = findpeaks(dPdt, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2-1));

% find peaks of inverted data
[~, MinIdx] = findpeaks(DataInv, 'MinPeakHeight',100, 'MinPeakDistance',length(dPdt)/(time_end*2-1));

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

%%  Find EDP: 0.2*dP/dt max
% Make sure the number of Minima = Maxima
if length(MinIdx) == length(pksT)
    
    % Pre allocate variables
    EDP = length(pksT);
    EDP_T = length(pksT);
    EDP_NT = length(pksT); % N for Negative (on the other side of pressure wave, negative derivative)
    EDP_N = length(pksT);
    BAD_First_WAVE = false;
    BAD_Last_WAVE = false;
    
    % scroll through all maxima
    for i = 1:length(pksT)
        
        EDi = pksT(i);
        while dPdt(EDi) > 0.2*pks(i)
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
                    BAD_Last_WAVE = true; % record flag that remarks a bad last wave
                end
            end
            
            if BAD_Last_WAVE == false
                % find which point is closest to the pressure value of EDP
                tempNeg_EDPs = EDP_Ni-3:1:EDP_Ni+3;

                % calculate the difference between the local neighborhood
                % (tempNeg_EDPs) and EDP
                EDP_Diffs = abs(EDP(i)-Pres(tempNeg_EDPs));

                % find minimum difference
                [~, tempInds] = min(EDP_Diffs);

                % assign Negative EDP values
                EDP_NT(i) = tempNeg_EDPs(tempInds); % keep in mind this is an indicie of the time vector
                EDP_N(i) = Pres(tempNeg_EDPs(tempInds)); % pressure in mmHg
            end
                        
        end
    end
    
    % ****
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
    
else 
    disp('SOMETHING BAD!!!!!!!!!!!!!!');
    disp('Number of minima does not match number of maxima!!!!!!!!');
    disp('Check the plots to see what happened');
    % --------------------------------------------
    %%% Maybe make some UI to remove bad waves???
    
    % Discuss with Hunter and Melanie
    % --------------------------------------------
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

t=time; pressure=Pres;
n=length(t);

% pre-allocate
A = zeros(n/2,1);
B = zeros(n/2,1);
Fourier = zeros(n,1);
error = zeros(n,1);

for i=1:n
    
    for k=1:(n/2) 
        A(k)=(2/n)*pressure(i)*cos((i*(2*pi)*t(i))/n); %A coefficients 
    end
    Asum=2/n+sum(A); %sum of A coeffs
    
    for k=1:(n/2)
        B(k)=(2/n)*pressure(i)*sin((i*(2*pi)*t(i))/n); % B coefficients
    end
    Bsum=sum(B); %sum of B coeffs
    
    Fourier(i)=(Asum*cos((i*(2*pi)*t(i))/n)) + (Bsum*sin((i*(2*pi)*t(i))/n)); 
    error(i) = (pressure(i)-Fourier(i))/pressure(i); %this is the difference between the raw value and the fourier interpolation
    % disp([num2str(t(i)),'      ',num2str(pressure(i)),'      ',num2str(Fourier(i)),'      ',num2str(error(i))]);
end

if FOUR_PRES == true
    %looks exactly the same as original data
    figure, hold on
    plot(t,pressure,'r',t,Fourier,'k');
    set(gca,'fontsize',14);
    title('Pressure Waveform','FontSize',20);
    legend('Raw Data','Fourier Series');
    xlabel('time (s)','FontSize',18);
    ylabel('RV Pressure (mmHg)','FontSize',18);
end

%-------------------------------------------------------------------------
% Recall dicsusion about interpolating points between points of isovolumic
% region in order to increase the number of points used for the sinusoid
% fitting

% Discuss with Hunter and Melanie
%-------------------------------------------------------------------------

%% Pes = 30 ms prior to dp/dt min %%%%%%%%%%%%%
dtmin_30=time(MinIdx)-0.03; %30 ms prior to dp/dt min
%this next section is a lot trickier...
%finding the corresponding pressure values at the time of dp/dt min

% pre - allocate
diff = zeros(5,length(dtmin_30));
mindiff = zeros(length(dtmin_30),1);

% Subtracting 0.03 seconds results in a timepoint that does not exist (as
% the sampling rate is 0.004 s or 4 ms)
% this loop extracts the closest time point to the 30 ms subtraction 
for i=1:length(dtmin_30)
    tempGuess = round(dtmin_30(i)/0.004); % convert time value to index
    tempRng = tempGuess-2:1:tempGuess+2;
    diff(:,i) = dtmin_30(i)-time(tempRng);
    [mindiff(i), tempInds] =min(abs(diff(:,i))~=0);
    P_esTimes(i) = time(tempRng(tempInds));
    P_es(i) = Pres(tempRng(tempInds));
end

if ES_PTS == true
    figure, hold on;
    plot(time,Pres,'b',P_esTimes,P_es,'ro',time(MinIdx), Pres(MinIdx), 'ko');
    set(gca,'fontsize',14);
    title('Pes 30 [ms] prior to dP/dt', 'FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressure [mmHg]','FontSize',18);
    legend('Pressure','Pes points', 'dP/dt Min');
    box on;
    grid on;
    hold off;
end
% hm, probably the most correct method
%We know that E_es should be a little less than the max pressure value on the curve.

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

totIsoTimePoints = [];
totIsoPresPoints = [];

% scroll through the number of rows (pressure waves) in the
% structures: isovoltime and isovol
for i = 1:length(EDP)
    
    WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];
    
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; %this
    % equation is from Naeiji et al, single beat method of VVC
    c2=[11, 210, 8, -1];
    %       sin_fun=@(P)(1/2*P(1)*(1-cos(timepeaks2(:,i)*P(3)+P(2)))+P(4))-peaks2(:,i); %this is the equation from Takeuchi et al calculating Pmax,
          % pt= 1/2 * PmAX ([1-cos(wt+C)) + EDP
    % c0=[800 2 13 11]; %initial guesses for Pmax, angular frequency, phase shift angle, and end diastolic pressure
    [c,resnorm,residual]=lsqnonlin(sin_fun2,c2); %least squares fitting
%     Psine_RV2=(c(1)+c(2)*sin(c(3)*CurrWaveTs+c(4))); %this is for
%     % the first equation
%     %     Psine_RV(:,i)=(1/2*c(1)*(1-cos(timepeaks2(:,i)*c(3)+c(2)))+c(4)); %this is taking the returned  c values and plotting it
%     r_square2=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2; %rsquared values. We dont really need this since
     % its not necessarily a tight fit, I just wanted to see the values.
    P_max2(i)=c(1)+2*c(2); %first equation pmax, A+2B
    %or the other definition of pmax
    % P_max(i)=c(1);
    c_tot2(i,:)=c; %getting all the c values in a matrix
%     clear Psine_RV2 r_square2;

    % store the time points and pressure points in one array for easy
    % plotting 
    totIsoTimePoints = [totIsoTimePoints; WaveTs];
    totIsoPresPoints = [totIsoPresPoints; WavePs];
end
% plotting the isovolumic points on top of the pressure curve, and the sine
% wave estimating isovolumic beat
if ISOVOL_REG == true
    
    figure; hold on;
    plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints,'ro');
    set(gca,'fontsize',14);
    title('Isovolumetric Sinusoid Region','FontSize',20);
    ylabel('Pressure [mmHg]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    grid on;
    box on;
    hold on;

    % Attain the sinusoid fit for all points (so Pmax can be visualized
    for i = 1:length(EDP)

        % obtain the range of time of each peak
        interval = time(isovoltime(i).PosIso(1,1)):0.004:time(isovoltime(i).NegIso(end,1));

        % plug into Naeiji equation that was just solved for 
        FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));
        
        plot(interval, FitSinePres, 'k--');
        hold on;
    end
    
    legend('Pressure','Isovolumic Points', 'Sinusoid Fit');
    hold off;
end

% -------------------------------------------------------------------------
% Insert Nick's QC check for GUI approving individual sinusoidal fitting

% Or ...  Insert GUI that allows user to interactively change the ICs of
% individual pressure waves???

% Discuss with Hunter and Melanie
% -------------------------------------------------------------------------     

% ----------------------------------------------------------
% store AVG_Pes and AVG_Pmax and VVCRinvert2 for analysis
% ----------------------------------------------------------

%OKAY! here are the final values and we can FINALLY calculate VVCR.

% average Pes
AVG_Pes=mean(P_es); 

% average P_max
AVG_Pmax=mean(P_max2); 

%this is the eqation from Uyen Truongs VVCR paper
VVCR_UT=AVG_Pes/(AVG_Pmax-AVG_Pes); 

%this is kendalls
VVCRinvert=(AVG_Pmax/AVG_Pes)-1;