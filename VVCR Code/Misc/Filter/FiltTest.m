%% Load data and plot unfiltered
[Pres, dPdt, ~, Pnam, Pmrn, file, ~, ~] ...
    =  loadp_10_10_16('C:\Users\hunterk\Desktop\VVCR Code\','HA002019.&03.txt',100);

fdpdt = figure;
plot(dPdt);
hold on;
%plottools('on','plotbrowser'); % Uncomment to enable selection window

%% Looking at different parameters of elliptic filter

for i = 1: 1: 10
    n = 10;
    Rp = 5;
    Rs = 250;
    TimeStep = 1/250;
    
    Wp = TimeStep*i*n;
    Wp = 0.04+i*0.01;
    
    [b,a] = ellip(n,Rp,Rs,Wp);
    dPdt_filt = filtfilt(b,a,dPdt);
    figure(fdpdt); plot(dPdt_filt);
    
    figure; freqz(b,a);
    % pause;
end

n = 10; Rp = 5; Rs = 100; Wp = 0.07;
n = 10; Rp = 5; Rs = 150; Wp = 0.08;
n = 10; Rp = 5; Rs = 200; Wp = 0.09;

% This was Lammers' approach to filter tissue stress-strain data, but it
% attenuates our peak dP/dt too much (qualitative assessment). Probably
% worked for stress-strain because there were no peaks - so it effectively
% removed all the HF noise.

%% Looking at Elliptic -vs- Butterworth
n = 20;
Rp = 5;
Rs = 100;
TimeStep = 1/250;

Wp = TimeStep*2*n

[b,a] = ellip(n,Rp,Rs,Wp);
dPdt_filt = filtfilt(b,a,dPdt);
figure(fdpdt); plot(dPdt_filt);
figure; freqz(b,a);

[b,a] = butter(n,Wp);
dPdt_filt = filtfilt(b,a,dPdt);
figure(fdpdt); plot(dPdt_filt);
figure; freqz(b,a);

% Rolloff slower in Butterworth but that's OK - zero bandpass attenuation,
% so no distortion of the signal we want to keep...

%% Looking at different parameters of Butterworth filter

for i = 0: 1: 8
    n = 10+5*i;
    TimeStep = 1/250;
    
    Wp = TimeStep*2*n;
    Wp = TimeStep*2*30;
    
    [b,a] = butter(n,Wp);
    dPdt_filt = filtfilt(b,a,dPdt);
    figure(fdpdt); plot(dPdt_filt);
    
    figure; freqz(b,a);
    % pause;
end
% n = 30, 35 seems to provide the best balance between noise loss and
% retention of peaks (for this sample set of RV waveforms).

% Holding Wp constant and increasing the order of the filter gives a more
% rapid rolloff after the passband. But it honestly makes little difference
% to the actual plot; filtered response looks the same - which suggests
% that the content being filtered out is at pretty high frequency anyway.
% There probably is no need to increase the order with the cutoff, but this
% is how MATLAB examples do it, so I'll stick with that...

%% Comparison of int(filt(dP/dt)) to true P

n = 30;
TimeStep = 1/250;
Wp = TimeStep*2*n;
disp(['Cutoff is ' num2str(1/(2*TimeStep)*Wp)]);

[b,a] = butter(n,Wp);
dPdt_filt = filtfilt(b,a,dPdt);
figure(fdpdt); plot(dPdt_filt);
%figure; freqz(b,a);

Pres_filt = cumtrapz(dPdt_filt)*TimeStep+Pres(1);

fpres = figure;
plot(Pres);
hold on;
plot(Pres_filt);

n = 35;
TimeStep = 1/250;
Wp = TimeStep*2*n;
disp(['Cutoff is ' num2str(1/(2*TimeStep)*Wp)]);

[b,a] = butter(n,Wp);
dPdt_filt = filtfilt(b,a,dPdt);
figure(fdpdt); plot(dPdt_filt);
axis([1520 1680 -Inf Inf]);    % "Smooth"
%axis([925 1075 -Inf Inf]);   % "Typical"
%axis([160 320 -Inf Inf]);    % "Noisy"

%figure; freqz(b,a);

Pres_filt = cumtrapz(dPdt_filt)*TimeStep+Pres(1);
figure(fpres);
plot(Pres_filt);
axis([1520 1680 -Inf Inf]);    % "Smooth"
%axis([925 1075 -Inf Inf]);
%axis([160 320 -Inf Inf]);

% Go with a 35Hz cutoff Butterworth. See how much this effects the analysis
% if the filtered data is actually used to get a fit... If it's too much,
% maybe use the filtered dP/dt for isovolumic relaxation but not filtered
% Pres for isovolumic contraction.