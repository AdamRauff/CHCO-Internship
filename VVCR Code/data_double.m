function [Ret] = data_double (Data_O, ivIdx)

% Copy all original data to returned structure.
Ret = Data_O;

% increase the number of data points to provide a more reliable fit of the
% isovolumic sinusoid / multiharmonic waveform.
mysz = length(Data_O.Pres)*2;
mystp = Ret.time_step/2;

Ret.Pres_D = interpft(Data_O.Pres, mysz);
Ret.dPdt_D = interpft(Data_O.dPdt, mysz);
Ret.Time_D = mystp:mystp:Data_O.Time(end);

%% Pes = 30 ms prior to dp/dt min; get that time and convert to an integer. 
dtmin1_30 = Data_O.Time(ivIdx.dPmin1)-0.03;
dtmin2_30 = Data_O.Time(ivIdx.dPmin2)-0.03;

Ret.Pes1Times = Ret.Time_D(uint16(dtmin1_30/mystp)); 
Ret.Pes1 = Ret.Pres_D(uint16(dtmin1_30/mystp));
Ret.Pes2Times = Ret.Time_D(uint16(dtmin2_30/mystp)); 
Ret.Pes2 = Ret.Pres_D(uint16(dtmin2_30/mystp));

