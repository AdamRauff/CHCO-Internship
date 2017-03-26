%% kendalls old kid  
%UT 12 Jocelyn Smith, HA002019.&03
clc
clear all
close all


[FileName,PathName]=uigetfile;

[Pres, dPdt, Rvals, nam, mrn, marks]=loadp(PathName,FileName,100);

time_end=0.004*size(Pres); time=0.004:0.004:time_end;


figure(),plot(time,Pres); %view the original data
figure(),plot(time,dPdt); %view the original data

 %FOURIER SERIES INTERPOLATION
 t=time; pressure=Pres;
for i=1:1:length(t)
    n=length(t);
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
%looks good
%  figure(),plot(t,pressure,'r',t,Fourier,'k'),title('Pressure Waveform'),legend('Raw Data','Fourier Series'),xlabel('time (s)'),ylabel('RV Pressure (mmHg)');

%BUILT IN METHOD
%THIS METHOD SUCKS.
% Perform fft, look at magnitude info (sqrt(a_i^2+b_i^2))
n = length(pressure); xf = fft(pressure); xM = abs(xf)/n;
xa = 0;
for i = 1 : 1 : n/2
ft = xf(i+1)/n;
xa = xa + real(ft)*cos(i*t) - imag(ft)*sin(i*t);
end
xa = 2*xa + xM(1);
%Plot this on top of the original dataset.
% figure(),plot(t,pressure,'bo','LineWidth',2,'MarkerSize',5); hold on;
% plot(t,pressure,'k',t,xa,'r','LineWidth',2);

%this section is to isolate the peaks to fit the sinusoidal function
%isolating the peaks in the pressure curve. I chose points manually per
%patient
  jj=1; a=1; aa=1;
for ii=1:length(t)-1
%     deriv(ii)=Fourier(ii+1)-Fourier(ii)/t(ii+1)-t(ii); %use the interpolation values for smoothness
%  if abs(deriv)>20
%      peaks1(jj)=Fourier(ii+1);
%      timepeaks1(jj)=t(ii+1);
%        jj=jj+1;
%  end
 %trying a new method, not derivatives but using EDP
  if Fourier(ii)>10 %mmHg, end diastolic pressure
     peaks2(a,aa)=Fourier(ii); %getting values for each peak in a column of m length with n # of peaks
     timepeaks2(a,aa)=t(ii); %getting values of time for each peak in a column of m length with n # of peaks
       a=a+1; %increase the row (continuing on down the column)
  end
 if Fourier(ii)<10 && Fourier(ii+1)>10 %this if statement is to increase onto the next peak when the pressure
     %curve rises above EDP again, i.e. the start of a new peak.
         aa=aa+1; %moves onto the next peak
         a=1;%This is to start the row back to one for the next peak
 end
end

% figure(),plot(timepeaks2,peaks2,'o'); %looks good

%SINUSOIDAL FITTING
for i=1:size(peaks2,2)
EDP=10;
sin_fun=@(P)(P(1)+P(2)*sin(P(3)*timepeaks2(:,i)+P(4)))-peaks2(:,i); %this
% equation is from Naeiji et al, single beat method of VVC
      c0=[1, 200, 8, -1]; %may change to have a better fitting
%       sin_fun=@(P)(1/2*P(1)*(1-cos(timepeaks2(:,i)*P(3)+P(2)))+P(4))-peaks2(:,i); %this is the equation from Takeuchi et al calculating Pmax,
      % pt= 1/2 * PmAX ([1-cos(wt+C) + EDP
% c0=[800 2 13 11]; %initial guesses for Pmax, angular frequency, phase shift angle, and end diastolic pressure
   [c,resnorm,residual]=lsqnonlin(sin_fun,c0); %least squares fitting
 Psine_RV(:,i)=(c(1)+c(2)*sin(c(3)*timepeaks2(:,i)+c(4))); %this is for
% the first equation
%     Psine_RV(:,i)=(1/2*c(1)*(1-cos(timepeaks2(:,i)*c(3)+c(2)))+c(4)); %this is taking the returned  c values and plotting it
 r_square(i)=1-resnorm/norm(Psine_RV(:,i)-mean(Psine_RV(:,i)))^2; %rsquared values. We dont really need this since
 % its not necessarily a tight fit, I just wanted to see the values.
      P_max(i)=c(1)+2*c(2); %first equation pmax, A+2B
%or the other definition of pmax
% P_max(i)=c(1);
c_tot(i,:)=c; %getting all the c values in a matrix
end
% figure(),plot(t,pressure,timepeaks2,Psine_RV,'*'),title('Sinusoidal fitting to peaks'), ylabel('Pressure (mmHg)'), xlabel('Time (s)');
%EDP = 11 mmHg
%take the 5 highest peaks
Pmax_sort=sort(P_max,'descend'); Pmax_top5=Pmax_sort(1:6);

%NOTES
% need to find End systolic pressure. Several definitions:
% 30 ms before dP/dt min. or intersection of sine curve with pressure
% curve. OR Max of the pressure curve. Or mPap, which is really inaccurate. hm...


%This is for determining END SYSTOLIC PRESSURE using the 30 ms prior to
%dp/dt min
for I=1:length(Fourier)-1
    deriv(I)=(Fourier(I+1)-Fourier(I))/(t(I+1)-t(I)); %this is dP/dt, using fourier values
end
timeplot=t(1:end-1);

deriv=deriv';
DataInv = 1.01*max(deriv) - deriv;
[Minima1,MinIdx] = findpeaks(DataInv);
  Minima = deriv(MinIdx); %finding the local minima in the derivative 
  rty=1;
  for bf=1:length(Minima)
      if Minima(bf)<-305 && Minima(bf)>-4000 %limiting minima to just the lowest minimums but eliminating the extreme errors
          trueMinima(rty)=Minima(bf); %put the true minimums (dp/dt min) in a vector
          rty=1+rty; %increase the vector
      end
  end
  
  q=1;
for y=1:length(timeplot)
    for yy=1:length(trueMinima) 
    if deriv(y) == trueMinima(yy) %FINDING  the corresponing time for the true minima (dp/dt min)
        timeminall(q)= timeplot(y); %put the corresponding time points in a vector
        q=q+1; %increase vector
    end
    end
end

l=1;
for yyy=1:length(Fourier)
    for yh=1:length(timeminall)
        if t(yyy)==timeminall(yh) 
            dpdtminpressureall(l)=Fourier(yyy); %finding the correpsonding pressure value for dp/dt min
            l=l+1;
        end
    end
end
q1=1;
for ass=1:length(dpdtminpressureall)
   if dpdtminpressureall(ass)>23
       dpdtminpressure(q1)=dpdtminpressureall(ass);
       timemin(q1)=timeminall(ass);
       q1=q1+1;
   end
end

%%%%%%%%%%%%%%%% NOW, finding dP/dt max, according to Naeji's method%%%%%%%%%%%%%%%
[pksmax,locsmax]= findpeaks(deriv,timeplot);
fff=1;
for ff=1:length(pksmax)
    if pksmax(ff)>386
        pkstruemax(fff)=pksmax(ff);
        locstruemax(fff)=locsmax(ff);
        fff=fff+1;
    end
end
  q=1;
for y=1:length(timeplot)
    for yy=1:length(pkstruemax) 
    if deriv(y) == pkstruemax(yy) %FINDING  the corresponing time for the true maximum in the derivative (dp/dt max)
        timemaxall(q)= timeplot(y); %put the corresponding time points in a vector
        q=q+1; %increase vector
    end
    end
end
l=1;
for yyy=1:length(Fourier)
    for yh=1:length(timemaxall)
        if t(yyy)==timemaxall(yh) %if a time point is equal to the time point at dP/dt max
            dpdtmaxpressureall(l)=Fourier(yyy); %finding the correpsonding pressure value for dp/dt max at that time point
            l=l+1;
        end
    end
end
%%%%% gotta clean up the dp/dtmax and get rid of double points 
ddd=1;
for vvv=1:length(dpdtmaxpressureall)
    if dpdtmaxpressureall(vvv)<21
        dpdtmaxpressure(ddd)=dpdtmaxpressureall(vvv);
        timemax(ddd)=timemaxall(vvv);
        ddd=ddd+1;
    end
end


        
figure(),plot(timeplot,deriv,'b',timeminall,trueMinima,'ro',locstruemax,pkstruemax,'go'),title('Derivative over time'),legend('Derivative','dP/dt min','dP/dt max');
figure(),plot(t,Fourier,'b',timemin,dpdtminpressure,'ro',timemax,dpdtmaxpressure,'go'),title('dP/dt min on the pressure curve'),legend('RV pressure','dP/dt min points','dP/dt max points');

%%%%%%%% 30 ms prior to dp/dt min %%%%%%%%%%%%%
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
 for ff=1:length(P_es)
     if P_es(ff)>25 %only concerned about high peaks and realistic Pes
         P_es_true(vb)=P_es(ff);
         dtmin_30_true(vb)=dtmin_30(ff);
         vb=vb+1;
     end
 end
 

% NOW PLOT IT TO SEEEEE

figure(),plot(t,Fourier,'b',dtmin_30_true,P_es_true,'ro'),xlabel('Time (s)'),ylabel('Pressure (mmHg)'),legend('Original Pressure curve','End Systolic Pressure points')
% hm, probably the most correct method
%We know that E_es should be a little less than the max pressure value on the curve.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%this is to fit the sinusoidal curve, P= a+b*sin(c*t+d) to the
%%%%%%isovolumetric regions of the pressure curve, i.e. from end diastolic
%%%%%%pressure to dP/dt max, and from dP/dt min to end diastolic pressure



%%%%%%first, getting EDP from 200 mmhg/s thresh-hold in the dP/dt curve. so
%%%%%%much analyzing....ugh
dd=1;
for nn=1:length(deriv)
   closeto200(nn)=((deriv(nn))-200)/200;
   if abs(closeto200(nn))<10^-1
       time200(dd)=timeplot(nn);
       dd=dd+1;
   end
end
closeto200=closeto200';
rr=1;
for dd=1:length(Fourier)-1
    for ss=1:length(time200)
    if timeplot(dd)==time200(ss)
    pres200(rr)=Fourier(dd);
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
    if Fourier(gg)>9 && Fourier(gg)<dpdtmaxpressure(hh) %if the Pressure value is above EDP and below dP/dt max
        isovol(bh,kk)=Fourier(gg); %put that pressure value in a vector
        isovoltime(bh,kk)=t(gg); %put the corresponding time in a vector
        bh=bh+1; %increase the vector
    end
         end
         if deriv2(gg)<0 %right side of the peak, descending limb
            if Fourier(gg)>9 && Fourier(gg)<dpdtminpressure(hh) %if the pressure value is below dp/dt min and above EDP
        isovol(bh,kk)=Fourier(gg); %add it to the vector, still the same column
        isovoltime(bh,kk)=t(gg); %add it to the time vector, still the same column as above
        bh=bh+1;
            end
         end
            if Fourier(gg)>9 && Fourier(gg+1)<9 %if the pressure value dips below EDP, move onto the next peak, i.e. move onto the next column
        kk=kk+1; bh=1;%restart the row back at 1 for a new column
            end
                if Fourier(gg)==dpdtminpressure(hh) %%when the pressure equals dpdt min, move onto the next dpdt min and max
         hh=hh+1;
          if hh>length(dpdtminpressure)|| hh>length(dpdtmaxpressure)
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
        

%OKAY! here are the final values and we can FINALLY calculate VVCR.
AVG_Pes=mean(P_es_true); %taking the average of the Pes
 AVG_Pmax=mean(Pmax_top5); AVG_Pmax2=mean(Pmax_top52);%taking the average of Pmax
 VVCR=AVG_Pes/(AVG_Pmax-AVG_Pes); %this is the eqation from Uyen Truongs VVCR paper
VVCRinvert=(AVG_Pmax/AVG_Pes)-1; %this is kendalls
 VVCRinvert2=(AVG_Pmax2/AVG_Pes)-1; %this is kendalls