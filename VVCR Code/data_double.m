function [Ret] = data_double (Data_O, ivIdx)

% TIME STEP 0.002 USED HERE, NEEDS TO BE GENERALIZED TO DATA TYPE (HUMAN/CALF)

% Copy all original data to returned structure.
Ret = Data_O;

% increase the number of data points to provide a tighter fit of the 
% sinusoid that is fitted to the iv1Presumic points
Ret.Pres_D = interpft(Data_O.Pres,length(Data_O.Pres)*2);
Ret.Time_D = 0.002:0.002:Data_O.Time(end);

%% Pes = 30 ms prior to dp/dt min; get that time and convert to an integer. 
dtmin_30 = Data_O.Time(ivIdx.dPmin)-0.03;

Ret.P_esTimes = Ret.Time_D(uint16(dtmin_30/0.002)); 
Ret.P_es = Ret.Pres_D(uint16(dtmin_30/0.002));

