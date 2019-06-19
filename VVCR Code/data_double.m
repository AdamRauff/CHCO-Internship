function [Ret] = data_double (Data_O, ivIdx)
% This function doubles the size of the data using a Fouier-transform based
% interpolation technique (nice job Adam! It's periodic data after all!).
% All fitting is performed on this data, which has a timestep of 2ms for
% clinical, 0.5ms for calf data.
% It also finds the Pes point at 30ms prior to (dP/dt)min (the "dog" method
% of finding Pes). The Vanderpool way (Pes3) is found previously in
% data_isoidx; it's merely taken into the doubled pointspace here.

% Copy all original data to returned structure.
Ret = Data_O;

% increase the number of data points to provide a more reliable fit of the
% isovolumic sinusoid / multiharmonic waveform.
mysz = length(Data_O.Pres)*2;
mystp = Ret.time_step/2;

Ret.Pres_D = interpft(Data_O.Pres, mysz);
Ret.dPdt_D = interpft(Data_O.dPdt, mysz);
Ret.dP2t_D = interpft(Data_O.dP2t, mysz);
Ret.Time_D = mystp:mystp:Data_O.Time(end);

%% Pes approximations.
% Pressure acceleration: Pes is at the maximum of d2P/dt2 ajust before
% (dP/dt)min. The finding technique code below is from data_isoseg (which
% actually is called for the rest of the times just after this routine). Note
% that these are VALUES, time and pressure. NOT INDICES. (RIGHT?)
mysz = length(ivIdx.PesP);
Ret.PesPTimes = zeros(mysz,1);
Ret.PesP = zeros(mysz,1);

for i = 1: 1: mysz
    % Find corresponding index of Pes in doubled data.
    Pes = find(round(Ret.Time_D,3) == round(Data_O.Time(ivIdx.PesP(i)),3));
    
    % Set values of time and pressure based on that index.
    Ret.PesPTimes(i) = Ret.Time_D(Pes);
    Ret.PesP(i) = Ret.Pres_D(Pes);
end
