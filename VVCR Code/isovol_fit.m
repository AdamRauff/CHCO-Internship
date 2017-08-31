function [RetVal] = isovol_fit (isovolPres, isovolTime, timeDoub, PresDoub, ICS)
%
% isovolPres      - Struct; Fitting values for contract and relax;
%                     isovolPres(i).PosIso - Pressure Values on contract
%                     isovolPres(i).NegIso - Pressure Values on relaxation
% isovolTime      - Struct; indices of contract and relax;
%                     isovolTime(i).PosIso - indices on contract
%                     isovolTime(i).NegIso - indices on relaxation
% timeDoub        - Vector; Time for entire dataset in doubled data
% PresDoub        - Vector; Pressure for entire dataset in doubled data
% ICS             - Struct or Vector; If called from VVCR_ (first call), is
%                   structure needed for individual-cycle ICs; if called from
%                   a GUI, contains contstant initial conditions for fit.

%% SINUSOIDAL FITTING
%%%%%%%finally, fitting sinusoid to isovolPresumetric regions%%%%%%%%%

opts1 = optimset ('Display', 'off');
nfits = length(isovolPres);

% Variables for main fit
c_tot2 = zeros(nfits,4);
P_max2 = zeros(nfits,1);
waveFit = zeros(nfits,1);
r_square2 = zeros(nfits,1);

% Ploting vectors of the fitting data for GUI_SINU_FIT  
totIsoTimePoints = [];
totIsoPresPoints = [];

% Variables for adding points to SINU_GUI plots (within Vanderpool method)
WHILE_LOOP_FLAG = zeros(nfits,1);
ADD_TPoints = []; 
ADD_PPoints = []; 

% scroll through the number of rows (pressure waves) in the
% structures: isovolTime and isovolPres
for i = 1:nfits
    
    WaveTs = [timeDoub(isovolTime(i).PosIso)'; timeDoub(isovolTime(i).NegIso)'];
    WavePs = [isovolPres(i).PosIso; isovolPres(i).NegIso];
    
    % this equation is from Naeiji et al, single beat method of VVC
    sin_fun2 = @(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 
    
    if ~isstruct (ICS)
        % ICs passed in from GUI
        c2 = ICS;
    else
        % deriving the initial values from the data
        % mean - average pressure value between dp/dt max and min (top of
        % curve)
        T1 = ICS.dPmaxIdx(i);
        T2 = ICS.dPminIdx(i);
        Mea = mean(double(ICS.Pres(T1:T2)));
    
        % Amplitude is twice the mean
        Amp = double(1.8*Mea);
    
        % keep in mind this means the initial conditions of every wave fit may
        % be slightly different, While values entered via GUI make ICs same for
        % all waves.
        c2 = [Mea, Amp, ICS.Freq, -0.5];
    end

    [c,resnorm,~]=lsqnonlin(sin_fun2,c2,[],[],opts1); %least squares fitting
    
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

    % store the time points and pressure points in one array for easy
    % plotting - first pass (call from VVCR_); otherwise, reconsitute these
    % arrays if needed just outside this loop.

    totIsoTimePoints = [totIsoTimePoints; WaveTs];
    totIsoPresPoints = [totIsoPresPoints; WavePs];
   
    % AR 6/5/17 -----------------------------------------------
    % adding points succesively to beginning of systole to make better fit
    % of sick patients with wide curves
    
    % obtain maximum pressure point on actual curve
    PresMax = max(PresDoub(isovolTime(i).PosIso(1,1):1:isovolTime(i).NegIso(end,1)));
    if r_square2(i) > 0.80 && P_max2(i) < PresMax
       
        % keep count of how many points added to systole side
        count = 0;

        temp_ADD_TPoints = [];
        temp_ADD_PPoints = [];

        while P_max2(i) < PresMax
            
            % add point to isovolTime(i).PosIso and corresponding isovolPres(i).PosIso
            isovolTime(i).PosIso = [(isovolTime(i).PosIso(1,1))-1; isovolTime(i).PosIso];
            isovolPres(i).PosIso = [PresDoub(isovolTime(i).PosIso(1,1)); isovolPres(i).PosIso];

            temp_ADD_TPoints = [ADD_TPoints; timeDoub(isovolTime(i).PosIso(1,1))];
            temp_ADD_PPoints = [ADD_PPoints; PresDoub(isovolTime(i).PosIso(1,1))];

            % update Wave(x)s variables
            WaveTs = [timeDoub(isovolTime(i).PosIso)'; timeDoub(isovolTime(i).NegIso)'];
            WavePs = [isovolPres(i).PosIso; isovolPres(i).NegIso];

            % mark flag that points are added
            WHILE_LOOP_FLAG(i) = true;

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
               WHILE_LOOP_FLAG(i) = false;
            else
               waveFit(i) = 0;
               WHILE_LOOP_FLAG(i) = true;
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
                disp('    isovol_fit: Added nine points on systolic side of curve, and Pmax');
                disp('        remains short of actual pressure');
                disp(['        Wave: ',num2str(i), 'is excluded']);
                break
            end
        end

        if WHILE_LOOP_FLAG(i) == true
            ADD_TPoints = [ADD_TPoints; temp_ADD_TPoints];
            ADD_PPoints = [ADD_PPoints; temp_ADD_PPoints];
        end
    end
    
    % ---------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 (last pressure waveform) for an example 
    
    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which can make a good fit (r^2 > 0.99).
    % --------------------------------------------------------------------
end

%% if iso points have been added, re-compose the totIsoPnts variables
if any(WHILE_LOOP_FLAG)

    totIsoTimePoints = [totIsoTimePoints; ADD_TPoints];
    totIsoPresPoints = [totIsoPresPoints; ADD_PPoints];

end

% Fill out return structure - passed to the GUIs or used to update the GUI
% global handles.

% These vars may be passed in but may be updated by Vanderpool section.
RetVal.ivTime_D = isovolTime;          % Fitting Points
RetVal.ivPres_D = isovolPres;

% This data is generated by the fit itself, so must be stored here.
RetVal.ivPlotTime  = totIsoTimePoints;  % Plotting Points
RetVal.ivPlotPres  = totIsoPresPoints;
RetVal.fit.BadCyc  = waveFit;           % Which waveforms had a bad fit
RetVal.fit.InitIC  = c2;                % First Intial conditions used
RetVal.fit.RCoef   = c_tot2;            % First regression constants
RetVal.fit.PIsoMax = P_max2;            % Pmax values obtained from fit

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(waveFit==1); % find indices of the bad waves
if ~isempty(indX)
    disp('    isovol_fit: The following waves did NOT have a good fit (will not be included)');
    disp(['        Wave(s): ', num2str(indX')]);
else
    disp('    isovol_fit: All waves seemed to fit well!');
end
