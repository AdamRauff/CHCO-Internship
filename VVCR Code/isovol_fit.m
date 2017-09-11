function [Ret1, Ret2, Ret3] = isovol_fit (ivSeg, Data, ICS)
%
% ivSeg  - Struct of all pres and time fitting values:
%            iv1Pres/iv1Time/iv2Pres/iv2Time 1st level structs; Time labels
%              actually contain indices, not real times.
%            PosIso/NegIso 2nd level structs; values for contract and relax
% Data   - Structure containing Time and Pressure.
% ICS    - Struct or Vector; If called from VVCR_ (first call), is struct
%          needed for individual-cycle ICs; if called from a GUI, contains
%          contstant initial conditions for fit.

% Set existing ivSeg to minimally insure continuity if nothing changes. 
Ret2 = ivSeg;

% Values for Kind fit; note that Pres actually contains P and dP/dt values.
% THESE WILL COME IN WITH ivSeg STRUCTURE.
%iv2Pres = ivSeg.iv2Pres;
%iv2Time = ivSeg.iv2Time; 

opts1 = optimoptions (@lsqnonlin);
opts1.Display = 'off';
opts1.MaxFunctionEvaluations = 2000;
opts1.MaxIterations = 1000;

% Variables for main fit
nfits = length(ivSeg.iv1Pres);
c_tot2 = zeros(nfits,4);
P_max2 = zeros(nfits,1);
waveFit = zeros(nfits,1);
r_square2 = zeros(nfits,1);

% Ploting vectors of the fitting data for GUI_SINU_FIT  
Ret3.ivPlotTime = [];
Ret3.ivPlotPres = [];

% Variables for adding points to SINU_GUI plots (within Vanderpool method)
VanderCyc = zeros(nfits,1);
ADD_TPoints = []; 
ADD_PPoints = []; 

% scroll through the number of rows (pressure waves) in the
% structures: ivSeg.iv1Time and ivSeg.iv1Pres
for i = 1:nfits
    
    WaveTs = [Data.Time_D(ivSeg.iv1Time(i).PosIso)'; ...
        Data.Time_D(ivSeg.iv1Time(i).NegIso)'];
    WavePs = [ivSeg.iv1Pres(i).PosIso; ivSeg.iv1Pres(i).NegIso];
    
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

    [c,resnorm,~] = lsqnonlin (sin_fun2,c2,[],[],opts1);
    
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

    Ret3.ivPlotTime = [Ret3.ivPlotTime; WaveTs];
    Ret3.ivPlotPres = [Ret3.ivPlotPres; WavePs];
   
    % AR 6/5/17 -----------------------------------------------
    % adding points succesively to beginning of systole to make better fit
    % of sick patients with wide curves
    
    % obtain maximum pressure point on actual curve
    PresMax = max(Data.Pres_D(ivSeg.iv1Time(i).PosIso(1,1):1: ...
        ivSeg.iv1Time(i).NegIso(end,1)));
    if r_square2(i) > 0.80 && P_max2(i) < PresMax
       
        % keep count of how many points added to systole side
        count = 0;

        temp_ADD_TPoints = [];
        temp_ADD_PPoints = [];

        while P_max2(i) < PresMax
            
            % add point to iv1Time(i).PosIso and iv1Pres(i).PosIso
            ivSeg.iv1Time(i).PosIso = ...
                [(ivSeg.iv1Time(i).PosIso(1,1))-1; ivSeg.iv1Time(i).PosIso];
            ivSeg.iv1Pres(i).PosIso = ...
                [Data.Pres_D(ivSeg.iv1Time(i).PosIso(1,1)); ...
                ivSeg.iv1Pres(i).PosIso];

            temp_ADD_TPoints = ...
                [ADD_TPoints; Data.Time_D(ivSeg.iv1Time(i).PosIso(1,1))];
            temp_ADD_PPoints = ...
                [ADD_PPoints; Data.Pres_D(ivSeg.iv1Time(i).PosIso(1,1))];

            % update Wave(x)s variables
            WaveTs = [Data.Time_D(ivSeg.iv1Time(i).PosIso)'; ...
                Data.Time_D(ivSeg.iv1Time(i).NegIso)'];
            WavePs = [ivSeg.iv1Pres(i).PosIso; ivSeg.iv1Pres(i).NegIso];

            % re-fit sinusiod
            % equation from Naeiji et al, single beat method of VVC
            sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

            %least squares fitting
            [c,resnorm,~]=lsqnonlin(sin_fun2,c2); 

            Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

            % r^2 value
            r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

            % if the fit of the wave was bad, mark that wave. Only mark that
            % points are added if we actually take the result.
            if r_square2(i) <0.90
               waveFit(i) = 1;
               VanderCyc(i) = 0;
            else
               waveFit(i) = 0;
               VanderCyc(i) = 1;
            end
            
            %getting all the c values in a matrix
            c_tot2(i,:)=c; 
    
            %first equation pmax, A+B
            P_max2(i)=c(1)+abs(c(2));
            
            % increment count to keep track of added points
            count = count +1;
            
            % Do not let program add more than 10 points
            if count >= 10 && (P_max2(i) < PresMax || waveFit(i) == 1)


                disp('    isovol_fit: Added nine points on systolic side of curve, and Pmax');
                disp('        remains short of actual pressure');
                disp(['        Wave: ',num2str(i), 'is excluded']);

                waveFit(i) = 1;
                VanderCyc(i) = 0;

                break
            end
        end

        if VanderCyc(i) == true
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
if any(VanderCyc)

    Ret3.ivPlotTime = [Ret3.ivPlotTime; ADD_TPoints];
    Ret3.ivPlotPres = [Ret3.ivPlotPres; ADD_PPoints];

    temp = 1:1:nfits;
    disp(['    isovol_fit: Vanderpool Points added on cycles ' ...
        num2str(temp(logical(VanderCyc)))]);

    % Update ivSeg; iv1Time & iv1Pres may be updated in Vanderpool section.
    Ret2 = ivSeg;

end

% Fill out return structure - passed to the GUIs or used to update the GUI
% global handles.

% Fit data, this should be stored in its own fit section.
Ret1.BadCyc  = waveFit;   % Which waveforms had a bad fit
Ret1.InitIC  = c2;        % First Intial conditions used
Ret1.RCoef   = c_tot2;    % First regression constants
Ret1.PIsoMax = P_max2;    % Pmax values obtained from fit
Ret1.VCyc    = VanderCyc; % Were points added to failing waveforms?

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
