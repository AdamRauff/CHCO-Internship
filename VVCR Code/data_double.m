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
% Original = 30ms prior to (dP/dt)min; get time and convert to an integer.
% Note that these are VALUES, time and pressure. NOT INDICES.
dtmin1_30 = Data_O.Time(ivIdx.dPmin1)-0.03;
dtmin2_30 = Data_O.Time(ivIdx.dPmin2)-0.03;

Ret.Pes1Times = Ret.Time_D(uint16(dtmin1_30/mystp));
Ret.Pes1 = Ret.Pres_D(uint16(dtmin1_30/mystp));
Ret.Pes2Times = Ret.Time_D(uint16(dtmin2_30/mystp));
Ret.Pes2 = Ret.Pres_D(uint16(dtmin2_30/mystp));

% Vanderpool, Pes is at the minimum of d2P/dt2 just before (dP/dt)min. The
% finding technique code below is from data_isoseg (which actually is
% called for the rest of the times just after this routine).

mysz3 = length(ivIdx.Ps3);
Ret.Pes3Times = zeros(mysz3,1);
Ret.Pes3 = zeros(mysz3,1);

for i = 1: 1: mysz3
    % Find corresponding index of Pes in doubled data.
    Pes = find(round(Ret.Time_D,3) == round(Data_O.Time(ivIdx.Pes3(i)),3));
    
    % Set values of time and pressure based on that index.
    Ret.Pes3Times(i) = Ret.Time_D(Pes);
    Ret.Pes3(i) = Ret.Pres_D(Pes);
end