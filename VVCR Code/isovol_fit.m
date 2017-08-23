function [ PeakStruct ] = isovol_fit ( EDP, EDP_Doub, EDP_NT_Doub, TotNumWaves, time_end, isovol, isovoltime, time, timeDoub, pksT, MinIdx, Pres, PresDoub )
%
% EDP         - Vector; End Disatolic Pressure Values (start of IsoCont)
% EDP_Doub    - Vector; Index of EDP in doubled data (end of IsoRelax)
% EDP_NT_Doub - Vector; Index of end of Isovolumetric Relaxation in doubled data
% TotNumWaves - Scalar; Total # of cycles
% time_end    - Scalar; Final time of dataset
% isovol      - Vector Structure; Fitting values for contract and relax;
%                   isovol(i).PosIso - Pressure Values on contract
%                   isovol(i).NegIso - Pressure Values on relaxation
% isovoltime  - Vector Structure; indices of contract and relax;
%                   isovoltime(i).PosIso - indices on contract
%                   isovoltime(i).NegIso - indices on relaxation
% time        - Vector; Time of entire dataset
% timeDoub    - Vector; Time of entire dataset in doubled data
% pksT        - Vector; Indicies of (dP/dt)_max
% MinIdx      - Vector; Indicies of (dP/dt)_min
% Pres        -
% PresDoub    -


%% SINUSOIDAL FITTING
%%%%%%%finally, fitting sinusoid to isovolumetric regions%%%%%%%%%

% pre - allocate 
c_tot2 = zeros(length(EDP),4);
P_max2 = zeros(length(EDP),1);
waveFit = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);

totIsoTimePoints = [];
totIsoPresPoints = [];

% freq is an initial condition that does not rquire individual wave
% calculation --> executed outside loop
% frequnecy is the conversion to angular frequency 2*pi/T
% multiplied by the number of waves found over the time period
Freq = double(((2*pi)*TotNumWaves)/(time_end));

% scroll through the number of rows (pressure waves) in the
% structures: isovoltime and isovol
for i = 1:length(EDP)
    
    WaveTs = [timeDoub(isovoltime(i).PosIso)'; timeDoub(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];
    
    % this equation is from Naeiji et al, single beat method of VVC
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 
    
    % deriving the initial values from the data
    % mean - average pressure value between dp/dt max and min (top of
    % curve)
    T1 = pksT(i);
    T2 = MinIdx(i);
    Mea = mean(double(Pres(T1:T2)));
    
    % Amplitude is twice the mean
    Amp = double(1.8*Mea);
    
    % keep in mind this means the initial conditions of every wave fit may
    % be slightly different, While values entered via GUI make ICs same for
    % all waves.
    c2=[Mea, Amp, Freq, -0.5];

    [c,resnorm,~]=lsqnonlin(sin_fun2,c2); %least squares fitting
    
    Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));
    
    % r^2 value
    r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;
    
    % if the fit of the wave was bad, mark that wave
    if r_square2(i) <0.90
       waveFit(i) = 1; 
    end
    
    %getting all the c values in a matrix
    c_tot2(i,:)=c; 
    
    %first equation pmax, A+B
    P_max2(i)=c(1)+abs(c(2)); 
   
    % AR 6/5/17 -----------------------------------------------
    % adding points succesively to beginning of systole to make better fit
    % of sick patients with wide curves
    
    % obtain maximum pressure point on actual curve
    PresMax = max(PresDoub(isovoltime(i).PosIso(1,1):1:isovoltime(i).NegIso(end,1)));
    if r_square2(i) > 0.80 && P_max2(i) < PresMax
       
        % keep count of how many points added to systole side
        count = 0;
        while P_max2(i) < PresMax
            
            % add point to isovoltime(i).PosIso and corresponding isovol(i).PosIso
            isovoltime(i).PosIso = [(isovoltime(i).PosIso(1,1))-1; isovoltime(i).PosIso];
            isovol(i).PosIso = [PresDoub(isovoltime(i).PosIso(1,1)); isovol(i).PosIso];

            % update Wave(x)s variables
            WaveTs = [timeDoub(isovoltime(i).PosIso)'; timeDoub(isovoltime(i).NegIso)'];
            WavePs = [isovol(i).PosIso; isovol(i).NegIso];

            % re-fit sinusiod
            % equation from Naeiji et al, single beat method of VVC
            sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

            %least squares fitting
            [c,resnorm,~]=lsqnonlin(sin_fun2,c2); 

            Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

            % r^2 value
            r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

            % if the fit of the wave was bad, mark that wave
            if r_square2(i) <0.90
               waveFit(i) = 1;
            else
                waveFit(i) = 0;
            end
            
            %getting all the c values in a matrix
            c_tot2(i,:)=c; 
    
            %first equation pmax, A+B
            P_max2(i)=c(1)+abs(c(2));
            
            % increment count to keep track of added points
            count = count +1;
            
            % Do not let program add more than 10 points
            if count >= 10 && (P_max2(i) < PresMax || waveFit(i) == 1)
                waveFit(i) = 1;
                disp('Added nine points on systolic side of curve, and Pmax remains short of actual pressure');
                disp(['Wave: ',num2str(i), 'is excluded']);
                break
            end
        end
    end
    
    % ---------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 (last pressure waveform) for an example 
    
    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which can make a good fit (r^2 > 0.99).
    % --------------------------------------------------------------------

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
