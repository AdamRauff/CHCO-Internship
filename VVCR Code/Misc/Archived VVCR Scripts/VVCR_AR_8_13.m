%% Kendalls old kid  
% UT 12 Jocelyn Smith, HA002019.&03
clc
close all
clear all

% FLAGS: boolean variables used as flags
ORIG_PRES = true; % original pressure vs. time plot
ORIG_DPDT = true; % original dP/dt vs. time plot
FOUR_PRES = false; % fourier interpolation of pressure curve
CALC_BUILT_IN_FOUR = false; % calculate fourier interpolation using built in method
FFT_PRES = false; % plot the built in fourier interploation
ISO_PEAKS = true; % plot the isolated peaks using EDP for separation
SINU_FIT = true; % plotting the sinusoidal fitting using ISO_PEAKS points
PRES_FILT = true; % plot the line of the pressure filter that used on the corresponding dP/dt minima
DPDT_MIN_MAX = true; % isolating minima and maxima of dP/dt
PRES_DPDT_MIN_MAX = true; % plotting dP/dt minima and maxima on pressure data
% temp_flag = true;

% obtain text file that contains the data
[FileName,PathName] = uigetfile('*.*');

% apply loadp function to load the data into matlab variables
[Pres, dPdt, Rvals, name, mrn, marks]=loadp(PathName,FileName,100);

% construct time array (units of 4 milliseconds from catheter machine)
time_end=0.004*size(Pres); time=0.004:0.004:time_end;

if ORIG_PRES == true
    % Plot original data
    figure(1),plot(time,Pres); hold on;
    set(gca,'fontsize',14);
    title('Original Pressure Vs. Time','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    box on
    grid on
    hold off;
end

if ORIG_DPDT == true
    % Plot derivative of original data
    figure(2),plot(time,dPdt); hold on;
    set(gca,'fontsize',14);
    title(' dP/dt Vs. Time','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    box on
    grid on
    hold off;
end

% UserInput = input('Press Enter to proceed');
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
    %looks good (as in looks exactly the same as original data)
    figure(),plot(t,pressure,'r',t,Fourier,'k'),title('Pressure Waveform')
    legend('Raw Data','Fourier Series'),xlabel('time (s)'),ylabel('RV Pressure (mmHg)');
end

% this method is not working properly as currently written
if CALC_BUILT_IN_FOUR == true
    %BUILT IN METHOD
    % Perform fft, look at magnitude info (sqrt(a_i^2+b_i^2))
    n = length(pressure); xf = fft(pressure); xM = abs(xf)/n; xM2 = abs(xf)/sqrt(n);
    xa = 0;
    for i = 1 : n/2
        ft = xf(i+1)/n;
        xa = xa + real(ft)*cos(i*t) - imag(ft)*sin(i*t);
    end
    xa = 2*xa + xM(1);
end

if FFT_PRES == true
    %Plot this on top of the original dataset.
    figure(), hold on;
    plot(t,pressure,'k',t,xa,'r','LineWidth',2);
    set(gca,'fontsize',14);
    title('Visulaizing built in Fourier interpolation','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    legend('Original','Fourier');
    box on;
    grid on;
    hold off;
end

%% Isolate Peaks

% Isolate the peaks to fit the sinusoidal function
% isolating the peaks in the pressure curve. I chose points manually per
% patient
 
rowInd=1; % row index 
colInd=1; % column index

%%%%%%%%
EDP=10; %mmHg, end diastolic pressure
%%%%%%%%

for i=1:length(t)-1
    
  % using EDP to isolate peaks
  if Fourier(i)>EDP 
     peaks2(rowInd,colInd)=Fourier(i); % getting values for each peak stacked in a column
     timepeaks2(rowInd,colInd)=t(i); % time values
       rowInd=rowInd+1; %increase row index
  end
  
  % when the pressure curve rises above EDP again, i.e. the start of a new peak.
  if Fourier(i)<10 && Fourier(i+1)>10 
         colInd=colInd+1; %moves column the next peak
         rowInd=1;% reset row index to 1
  end
  
end

% pre - allocate
PressPeakRng = zeros(2,size(timepeaks2,2));

% Obtain the time range of each peak
for i = 1:size(timepeaks2,2)
    
    PressPeakRng(1,i) = timepeaks2(1,i);
    
    if PressPeakRng(1,i) ~= 0
        
        % look for last nonzero entry
        for j = 1:size(timepeaks2,1)
            
            if timepeaks2(j,i) == 0 && timepeaks2(j-1,i) ~= 0
               PressPeakRng(2,i) = timepeaks2(j-1,i);
            elseif timepeaks2(j,i) ~=0 && timepeaks2(end,i) == timepeaks2(j,i)
               PressPeakRng(2,i) = timepeaks2(j,i);
            end
            
        end
    end
end

% plot the isolated peaks on top of the pressure data
if ISO_PEAKS == true
    figure; hold on;
    plot(t, pressure,'k', timepeaks2,peaks2,'o'); %looks good
    set(gca,'fontsize',14);
    title('Isolating Individual Peaks','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressure [mmHg]','FontSize',18);
    box on;
    grid on;
    hold off;
end

% UserInput = input('Press Enter to proceed');
%% SINUSOIDAL FITTING

% This method isn't primarily used. Ordinarily use the isovolumetric region
% for sinusoidal fitting.
% note this fitting is taking place with all points of each pressure
% waveform, rather than the isovolumetric regions

% pre -allocate variables
P_max = zeros(size(peaks2,2),1);
r_square = zeros(size(peaks2,2),1);

% Iterate through number of peaks
for i=1:size(peaks2,2)
    % Naeiji et al, single beat method of VVC
    sin_fun=@(P)(P(1)+P(2)*sin(P(3)*timepeaks2(:,i)+P(4)))-peaks2(:,i);
    
    % intial conditions
    c0=[1, 200, 8, -1]; %may change to have a better fitting
    
    % sin_fun=@(P)(1/2*P(1)*(1-cos(timepeaks2(:,i)*P(3)+P(2)))+P(4))-peaks2(:,i); %this is the equation from Takeuchi et al calculating Pmax,
    % pt= 1/2 * PmAX ([1-cos(wt+C) + EDP
    % c0=[800 2 13 11]; %initial guesses for Pmax, angular frequency, phase shift angle, and end diastolic pressure
    
    %least squares fitting
    [c,resnorm,residual]=lsqnonlin(sin_fun,c0); 
    
    % First equation
    Psine_RV(:,i)=(c(1)+c(2)*sin(c(3)*timepeaks2(:,i)+c(4))); 
    %     Psine_RV(:,i)=(1/2*c(1)*(1-cos(timepeaks2(:,i)*c(3)+c(2)))+c(4)); %this is taking the returned  c values and plotting it
    
    r_square(i)=1-resnorm/norm(Psine_RV(:,i)-mean(Psine_RV(:,i)))^2; %rsquared values. We dont really need this since
    % its not necessarily a tight fit, I just wanted to see the values.
    P_max(i)=c(1)+2*c(2); %first equation pmax, A+2B
    %or the other definition of pmax
    % P_max(i)=c(1);
    c_tot(i,:)=c; %getting all the c values in a matrix
end

if SINU_FIT == true
    % plot the sinusoidal fitting to each pressure curve using all of the
    % waveform points
    figure; hold on;
    plot(t,pressure,'k',timepeaks2,Psine_RV,'*');
    set(gca,'fontsize',14);
    title('Sinusoidal Fitting to Peaks','FontSize',20);
    ylabel('Pressure [mmHg]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    legend('Original Pressure','Best-fit Sinusoids');
    box on;
    grid on;
    hold off;
end

% UserInput = input('Press Enter to proceed');

%EDP = 11 mmHg
%take the 5 highest peaks
Pmax_sort=sort(P_max,'descend'); 
Pmax_top5=Pmax_sort(1:6);

%% FINDING dP/dt Minima

% need to find End systolic pressure (Pes). Several definitions:
% 30 ms before dP/dt min. or intersection of sine curve with pressure
% curve. OR Max of the pressure curve. Or mPap, which is really inaccurate. hm...

%This is for determining END SYSTOLIC PRESSURE using the 30 ms prior to
%dp/dt min

% ----------------------------------------------------------
% WHY NOT USE dPdt RETURNED FROM loaddp FUNCTION EARLIER ???

% answer: because the original data may often contain much noise.
% The newly calculated dP/dt (deriv) is calculated using the fourier
% interpolation values, which should give a smoother result. 
% -----------------------------------------------------------

% pre - allocate variable
deriv = zeros(length(Fourier)-1,1);

% This is a derivative that is taken from the fourier values
for i=1:length(Fourier)-1
    
    % dP/dt using fourier values
    deriv(i,1)=(Fourier(i+1)-Fourier(i))/(t(i+1)-t(i));  
end
timeplot=t(1:end-1);

% deriv=deriv';
% -----------------------------------
% WHY NOT JUST:
% DataInv = (-1)*deriv
% 
% answer: This is a way to flip the data and ensure all values are positive
% -----------------------------------

DataInv = 1.01*max(deriv) - deriv;

% use find peaks function to find the peaks of inverted data
[Minima1,MinIdx] = findpeaks(DataInv);

%finding the local minima values in deriv (dP/dt)
Minima = deriv(MinIdx); 

% scroll through minima and get rid of extreme errors
% May be subject to change per person

%%%%%%%%%%%%%%%%%%%%%%
dPdt_thresh = -305; % This is the minimum threshold value in order to identify a minima in dP/dt
%%%%%%%%%%%%%%%%%%%%%%

Index=1;
for i=1:length(Minima)
  if Minima(i)<dPdt_thresh && Minima(i)>-4000 %limiting minima to just the lowest minimums but eliminating the extreme errors
      trueMinima(Index)= Minima(i); %put the true minimums (dp/dt min) in a vector
      timeminall(Index) = timeplot(MinIdx(i)); %FINDING  the corresponing time for the true minima (dp/dt min)
      Index=1+Index; %increase the vector
  end
end

ind = 1;

% filter minima based on local range filtering
for i = 1:size(PressPeakRng,2)
    if PressPeakRng(1,i) ~= 0 && PressPeakRng(2,i) ~= 0
        % find all minima within the ith range (time)
        rngTime = find(timeminall > PressPeakRng(1,i) & timeminall < PressPeakRng(2,i));
        
        % find the dP/dt values (mmHg/s)
        tempMinima = trueMinima(rngTime);
        
        % store the value (mmHg/s) and time in another vector
        if length(tempMinima) > 1
            [trueMinima2(ind), tempInd] = min(tempMinima);
            timeminall2(ind) = timeminall(rngTime(tempInd));
            ind = ind + 1;
        elseif length(tempMinima) == 1
            trueMinima2(ind) = tempMinima;
            timeminall2(ind) = timeminall(rngTime);
            ind = ind + 1;
        end
    end
end

% pre - allocate
dpdtminpressureall = zeros(length(timeminall2),1);

% finding the correpsonding pressure value for dp/dt mins
Index=1;
for i=1:length(Fourier)
    for j=1:length(timeminall2)
        if t(i)==timeminall2(j) 
            dpdtminpressureall(Index)=Fourier(i);  
            Index=Index+1;
        end
    end
end

% filtering the pressure values that are greater than 23 mmHg.
% this value may be subject to change on a patient to patient basis
%------------------------------
% great choice of iterator
% -----------------------------
% pres_filter_dpdtMin = 23;
% 
% ind=1;
% for i=1:length(dpdtminpressureall)
%    if dpdtminpressureall(i) > pres_filter_dpdtMin
%        dpdtminpressureall(ind)=dpdtminpressureall(i);
%        timemin(ind)=timeminall(i);
%        ind=ind+1;
%    end
% end

if PRES_FILT == true
    % Visualize the pressure filter employed above
    figure; hold on
    plot(time, Pres', 'b', timeminall2, dpdtminpressureall, 'ko');
%     plot([time(1), time(end-1)], [pres_filter_dpdtMin, pres_filter_dpdtMin],'r--');
    set(gca,'fontsize',14);
    title('Finding dP/dt Minima','FontSize',20);
    ylabel('Pressure [mmHg]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    legend('Original Pressure','All dP/dt Minima');
    box on;
    grid on;
    hold off;
end

% UserInput = input('Press Enter to proceed');

%% Finding dP/dt Max
% finding dP/dt max, according to Naeji's method%%%%%%%%%%%%%%%

[pksmax,locsmax]= findpeaks(deriv,timeplot);
ind=1;

% may be subject to change on patient basis
PksMaxThresh = (-1)*dPdt_thresh;

% filter pks to be greater than a threshold
for i=1:length(pksmax)
    if pksmax(i) > PksMaxThresh
        pkstruemax(ind)=pksmax(i);
        locstruemax(ind)=locsmax(i);
        ind=ind+1;
    end
end

ind = 1;
% filter minima based on local range filtering
for i = 1:size(PressPeakRng,2)
    if PressPeakRng(1,i) ~= 0 && PressPeakRng(2,i) ~= 0
        % find all minima within the ith range (time)
        rngTime = find(locstruemax > PressPeakRng(1,i) & locstruemax < PressPeakRng(2,i));
        
        % find the dP/dt values (mmHg/s)
        tempMaxima = pkstruemax(rngTime);
        
        % store the value (mmHg/s) and time in another vector
        if length(tempMaxima) > 1
            [pkstruemax2(ind), tempInd] = max(tempMaxima);
            locstruemax2(ind) = locstruemax(rngTime(tempInd));
            ind = ind + 1;
        elseif length(tempMaxima) == 1
            pkstruemax2(ind) = tempMaxima;
            locstruemax2(ind) = locstruemax(rngTime);
            ind = ind + 1;
        end
    end
end

% ind=1;
% for i=1:length(timeplot)
%     for j=1:length(pkstruemax) 
%         if deriv(i) == pkstruemax(j) %FINDING  the corresponing time for the true maximum in the derivative (dp/dt max)
%             timemaxall(ind)= timeplot(i); %put the corresponding time points in a vector
%             ind=ind+1; %increase vector
%         end
%     end
% end

% find correpsonding pressure value for dp/dt max time point
ind=1;
for i=1:length(Fourier)
    for j=1:length(locstruemax2)
        if t(i)==locstruemax2(j) %if a time point is equal to the time point at dP/dt max
            dpdtmaxpressureall(ind)=Fourier(i); 
            ind=ind+1;
        end
    end
end

% clean up the dp/dtmax and get rid of double points
% may be subject to change on patient basis
% dPdtMaxPresThresh = 21;
% 
% ind=1;
% for i=1:length(dpdtmaxpressureall)
%     if dpdtmaxpressureall(i) < dPdtMaxPresThresh
%         dpdtmaxpressureall(ind) = dpdtmaxpressureall(i);
%         timemax(ind)=locstruemax(i);
%         ind=ind+1;
%     end
% end

% plotting dP/dt and circled minima and maxima
if DPDT_MIN_MAX == true
    figure, hold on;
    plot(timeplot,deriv,'b',timeminall2,trueMinima2,'ro',locstruemax2,pkstruemax2,'go');
%     % plot pressure peak neighborhoods (vertical lines)
%     for i = 1 : size(PressPeakRng,2)
%         if PressPeakRng(1,i) ~= 0 && PressPeakRng(2,i) ~= 0
%             plot([PressPeakRng(1,i) PressPeakRng(1,i)],[-400 400], 'k');
%             plot([PressPeakRng(2,i) PressPeakRng(2,i)],[-400 400], 'k--');
%         end
%     end
    set(gca,'fontsize',14);
    title('dP/dt Minima and Maxima','FontSize',20);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    legend('Derivative','dP/dt min','dP/dt max');
    box on;
    grid on;
    hold off;
end

% plotting dP/dt Maxima and minima on pressure data
if PRES_DPDT_MIN_MAX == true
    figure, hold on;
    plot(t,Fourier,'b',timeminall2,dpdtminpressureall,'ro',locstruemax2,dpdtmaxpressureall,'go');
    set(gca,'fontsize',14);
    title('dP/dt min on the pressure curve','FontSize',20);
    ylabel('Pressure [mmHg]','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    legend('RV pressure','dP/dt Min','dP/dt Max');
    box on;
    grid on;
    hold off;
end

%% 30 ms prior to dp/dt min %%%%%%%%%%%%%
dtmin_30=timemin-0.03; %30 ms prior to dp/dt min
%this next section is a lot trickier...
%finding the corresponding pressure values at the time of dp/dt min

for K=1:length(dtmin_30)
    for H=1:length(Fourier)
        diff(H,K)=(dtmin_30(K)-t(H))/dtmin_30(K); %this is the difference
        %between the 30 ms prior time from the time to determine the
        %minumum. Basically im finding where dtmin_30 matches the time in
        %the time vector
    end
    mindiff(K)=min(abs(diff(:,K))); %the minimum is where they match the best, ie
    %very small difference!
end
 
for K=1:length(dtmin_30)
    for H=1:length(Fourier)
    if abs(diff(H,K)) == mindiff(K) %when the difference value matches the minimum, I need
        %that place in time and location within the vector. I need its H
        %location
        P_es(K)=Fourier(H); % Put those corresponding pressure values 
        %in a vector!
    end
    end
end

vb=1;
 for i=1:length(P_es)
     if P_es(i)>25 %only concerned about high peaks and realistic Pes
         P_es_true(vb)=P_es(i);
         dtmin_30_true(vb)=dtmin_30(i);
         vb=vb+1;
     end
 end
 
% NOW PLOT IT TO SEEEEE

figure(),plot(t,Fourier,'b',dtmin_30_true,P_es_true,'ro'),xlabel('Time (s)'),ylabel('Pressure (mmHg)'),legend('Original Pressure curve','End Systolic Pressure points')
% hm, probably the most correct method
%We know that E_es should be a little less than the max pressure value on the curve.

%%%%%%this is to fit the sinusoidal curve, P= a+b*sin(c*t+d) to the
%%%%%%isovolumetric regions of the pressure curve, i.e. from end diastolic
%%%%%%pressure to dP/dt max, and from dP/dt min to end diastolic pressure

%%%%%%first, getting EDP from 200 mmhg/s thresh-hold in the dP/dt curve. so
%%%%%%much analyzing....ugh
ind=1;
for i=1:length(deriv)
   closeto200(i)=((deriv(i))-200)/200;
   if abs(closeto200(i))<10^-1
       time200(ind)=timeplot(i);
       ind=ind+1;
   end
end
closeto200=closeto200';

rr=1;

for i=1:length(Fourier)-1
    for j=1:length(time200)
        if timeplot(i)==time200(j)
            pres200(rr)=Fourier(i);
            rr=rr+1;
        end
    end
end

% figure(),plot(t,Fourier,time200,pres200,'o'),title('points at 200 mmhg/s'),legend('Data','time points at 200mmHg/s');
s12=200*ones(1,length(timeplot));
figure(),plot(t,Fourier,timeplot,deriv,'r',timeplot,s12,'g',time200,pres200,'o'),title('points at 200 mmhg/s'),legend('Data','time points at 200mmHg/s');


%%%%Next, another fourier interpolation to get smoother derivative value
n1 = length(pressure); xf = fft(pressure); xM = abs(xf)/n1;
M = floor((n1+1)/100); %here, I divide by 100 for increased smoothness
a0 = xf(1)/n1;
an = 2*real(xf(2:M))/n1;
a6 = xf(M+1)/n1;
bn = -2*imag(xf(2:M))/n1;
qt = 1:length(an); newt1=0:(t(end)/n1):t(end);
Fourier3 = a0 + an'*cos(2*pi*qt'*newt1/t(end)) ...
       + bn'*sin(2*pi*qt'*newt1/t(end)) ...
       + a6*cos(2*pi*6*newt1/t(end));
%  figure(),plot(t,Fourier,'b',newt1,Fourier3,'r');
%%%%%% finding the derivative of this new fourier interpolation
for I=1:length(Fourier3)-1
    deriv2(I)=(Fourier3(I+1)-Fourier3(I))/(newt1(I+1)-newt1(I)); %this is dP/dt, using fourier3 values
end
timeplot2=newt1(1:end-1);
 figure(),plot(timeplot2,deriv2);
%%%%%% fourier interp of the derivative

[pks,locs]=findpeaks(Fourier,t);
cc=1;
for v=1:length(pks)
    if pks(v)>70
        highpks(cc)=pks(v);
        locshigh(cc)=locs(v);
        cc=cc+1;
    end
end
% figure(),plot(t,Fourier,'b',locshigh,highpks,'o');
%%%%%%%%%THIS is to isolate the true isovolumetric region%%%%%%
bh=1; kk=1; hh=1;

for gg=1:length(Fourier)-1 %going through iterations of pressure first for every value of dtdtmaxpressure
    if deriv2(gg)>0 %left side of the peak, ascending limb
        if Fourier(gg)>9 && Fourier(gg)<dpdtmaxpressureall(hh) %if the Pressure value is above EDP and below dP/dt max
            isovol(bh,kk)=Fourier(gg); %put that pressure value in a vector
            isovoltime(bh,kk)=t(gg); %put the corresponding time in a vector
            bh=bh+1; %increase the vector
        end
    end
    if deriv2(gg)<0 %right side of the peak, descending limb
        if Fourier(gg)>9 && Fourier(gg)<dpdtminpressureall(hh) %if the pressure value is below dp/dt min and above EDP
            isovol(bh,kk)=Fourier(gg); %add it to the vector, still the same column
            isovoltime(bh,kk)=t(gg); %add it to the time vector, still the same column as above
            bh=bh+1;
        end
    end
    if Fourier(gg)>9 && Fourier(gg+1)<9 %if the pressure value dips below EDP, move onto the next peak, i.e. move onto the next column
        kk=kk+1; bh=1;%restart the row back at 1 for a new column
    end
    if Fourier(gg)==dpdtminpressureall(hh) %%when the pressure equals dpdt min, move onto the next dpdt min and max
        hh=hh+1;
        if hh>length(dpdtminpressureall)|| hh>length(dpdtmaxpressureall)
            break
        end
    end
end
 
 %%%%%%%need to clean up the data, because of the weird little peaks in the
 %%%%%%%ascending limbs and begining of diastole
  zz=1;
 for vv=1:size(isovol,2)
    for xx=1:size(isovol,1)-1
        if isovol(xx,vv)>0 && isovol(xx+1,vv)==0 || isovol(xx,vv)>0 && isovol(xx+1,vv)==isovol(end,vv)
             if xx>5
                 correct_isovol(:,zz)=isovol(:,vv);
                 correct_isovolt(:,zz)=isovoltime(:,vv);
                 zz=zz+1;
             end
        end
    end
 end
 
 %%%%%%%finally, fitting sinusoid to isovolumetric regions%%%%%%%%%
 %SINUSOIDAL FITTING
for i=1:size(correct_isovol,2)
EDP=10;
sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*correct_isovolt(:,i)+P(4)))-correct_isovol(:,i); %this
% equation is from Naeiji et al, single beat method of VVC
      c2=[9, 200, 8, -1];
%       sin_fun=@(P)(1/2*P(1)*(1-cos(timepeaks2(:,i)*P(3)+P(2)))+P(4))-peaks2(:,i); %this is the equation from Takeuchi et al calculating Pmax,
      % pt= 1/2 * PmAX ([1-cos(wt+C)) + EDP
% c0=[800 2 13 11]; %initial guesses for Pmax, angular frequency, phase shift angle, and end diastolic pressure
   [c,resnorm,residual]=lsqnonlin(sin_fun2,c2); %least squares fitting
 Psine_RV2(:,i)=(c(1)+c(2)*sin(c(3)*correct_isovolt(:,i)+c(4))); %this is for
% the first equation
%     Psine_RV(:,i)=(1/2*c(1)*(1-cos(timepeaks2(:,i)*c(3)+c(2)))+c(4)); %this is taking the returned  c values and plotting it
 r_square2(i)=1-resnorm/norm(Psine_RV2(:,i)-mean(Psine_RV2(:,i)))^2; %rsquared values. We dont really need this since
 % its not necessarily a tight fit, I just wanted to see the values.
      P_max2(i)=c(1)+2*c(2); %first equation pmax, A+2B
%or the other definition of pmax
% P_max(i)=c(1);
c_tot2(i,:)=c; %getting all the c values in a matrix
end
figure(),plot(t,pressure,correct_isovolt,Psine_RV2,'*'),title('Sinusoidal fitting to peaks'), ylabel('Pressure (mmHg)'), xlabel('Time (s)');
%EDP = 11 mmHg
%take the 5 highest peaks
Pmax_sort2=sort(P_max2,'descend'); Pmax_top52=Pmax_sort2(4:6);
 
figure(),plot(t,Fourier,'b',isovoltime,isovol,'o'),title('Original isovolumetric region');
figure(),plot(t,Fourier,'b',correct_isovolt,correct_isovol,'o'),title('Corrected isovolumetric region');       

% ----------------------------
% store AVG_Pes and AVG_Pmax
% and VVCRinvert2
% ----------------------------

%OKAY! here are the final values and we can FINALLY calculate VVCR.
AVG_Pes=mean(P_es_true); %taking the average of the Pes
AVG_Pmax=mean(Pmax_top5); AVG_Pmax2=mean(Pmax_top52);%taking the average of Pmax
VVCR_UT=AVG_Pes/(AVG_Pmax-AVG_Pes); %this is the eqation from Uyen Truongs VVCR paper
VVCRinvert=(AVG_Pmax/AVG_Pes)-1; %this is kendalls
VVCRinvert2=(AVG_Pmax2/AVG_Pes)-1; %this is kendalls