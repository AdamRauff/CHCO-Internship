function [Ret1, Ret2, Ret3] = fit_takeuchi (ivSeg, Data, ICS)
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

opts1 = optimoptions (@lsqnonlin);
opts1.Display = 'off';
opts1.MaxFunctionEvaluations = 2000;
opts1.MaxIterations = 1000;
opts1.FiniteDifferenceType = 'central';

% Variables for main fit, all returned in Ret1
nfits = length(ivSeg.iv1Pres);

Ret1.Rsq = zeros(nfits,1);     % Goodness of fit coefficients
Ret1.RCoef = zeros(nfits,4);   % Fit regression constants
Ret1.BadCyc = zeros(nfits,1);  % Which waveforms had a bad fit
Ret1.PIsoMax = zeros(nfits,1); % Pmax,iso values obtained from fit
if isstruct(ICS)
    Ret1.CycICs = zeros(nfits,4); % Saved cycle-specific ICS
end

% Variables for adding points to SINU_GUI plots (within Vanderpool method)
Ret1.VCyc = zeros(nfits,1);
ADD_TPoints = []; 
ADD_PPoints = []; 

% Ploting vectors of the fitting data for GUI_SINU_FIT  
Ret3.ivPlotTime = [];
Ret3.ivPlotPres = [];

% scroll through the number of rows (pressure waves) in the
% structures: ivSeg.iv1Time and ivSeg.iv1Pres
for i = 1:nfits
    
    % Obtain vectors of fitting time, pressure; then normalize time to start
    % at zero (to obtain consistent phase).
    WaveTs = [Data.Time_D(ivSeg.iv1Time(i).PosIso)'; ...
        Data.Time_D(ivSeg.iv1Time(i).NegIso)'];
    WavePs = [ivSeg.iv1Pres(i).PosIso; ivSeg.iv1Pres(i).NegIso];
    
    WaveTsNorm = WaveTs-WaveTs(1);

    % this equation is from Naeiji et al, single beat method of VVC
    sin_fun2 = @(P)(P(1)+P(2)*sin(P(3)*WaveTsNorm+P(4)))-WavePs; 
    
    if ~isstruct (ICS)
        % ICs passed in from GUI
        c2 = ICS;

        lb = [  0.0   0.0  0.5*ICS(3) -2*pi/3];
        ub = [500.0 500.0    2*ICS(3)   -pi/3];
    else
        % Deriving the initial values from the data. Freq is an average
        % found of the function. Mean is the specific average pressure
        % value between dp/dt max and min (top of curve) specific to this
        % cycle.
        T1 = ICS.dPmaxIdx(i);
        T2 = ICS.dPminIdx(i);
        Mea = mean(double(ICS.Pres(T1:T2)));
    
        % Amplitude is about twice the mean
        Amp = double(1.8*Mea);
    
        % keep in mind this means the initial conditions of every wave fit
        % may be slightly different, While values entered via GUI make ICs
        % same for all waves.
        c2 = [Mea, Amp, 1.5*ICS.Freq, -pi/2];
        Ret1.CycICs(i,:)= c2; % Saved cycle specific ICs

        lb = [  0.0   0.0    ICS.Freq -2*pi/3];
        ub = [500.0 500.0  2*ICS.Freq   -pi/3];
    end

    [c,SSE,~] = lsqnonlin (sin_fun2,c2,lb,ub,opts1);

    % r^2 value; if the fit was bad, mark that wave.
    SSTO = norm(WavePs-mean(WavePs))^2;
    Ret1.Rsq(i) = 1-SSE/SSTO;
    
    if Ret1.Rsq(i) <0.90
       Ret1.BadCyc(i) = 1; 
    end

    if any( (c-lb) < 0 ) || any ( (ub-c) < 0 )
       disp(['    fit_takeuchi: fit bounds violated on cycle ' ...
        num2str(i, '%02i')]);
       Ret1.BadCyc(i) = 1; 
    end
    
    % Store all coefficients, Pmax for return
    Ret1.RCoef(i,:) = c;
    Ret1.PIsoMax(i) = c(1)+abs(c(2));

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
    if Ret1.Rsq(i) > 0.80 && Ret1.PIsoMax(i) < PresMax
       
        % keep count of how many points added to systole side
        count = 0;

        temp_ADD_TPoints = [];
        temp_ADD_PPoints = [];

        while Ret1.PIsoMax(i) < PresMax
            
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

            WaveTsNorm = WaveTs-WaveTs(1);

            % re-fit sinusiod
            % equation from Naeiji et al, single beat method of VVC
            sin_fun2 = @(P)(P(1)+P(2)*sin(P(3)*WaveTsNorm+P(4)))-WavePs; 

            %least squares fitting
            [c,SSE,~] = lsqnonlin (sin_fun2,c2,lb,ub,opts1);

            % r^2 value; if the fit was bad, mark that wave.
            SSTO = norm(WavePs-mean(WavePs))^2;
            Ret1.Rsq(i) = 1-SSE/SSTO;

            % Only mark that points are added if we actually take the result.
            if Ret1.Rsq(i) <0.90
               Ret1.BadCyc(i) = 1;
               Ret1.VCyc(i) = 0;
            else
               Ret1.BadCyc(i) = 0;
               Ret1.VCyc(i) = 1;
            end
            
            % Store all coefficients, Pmax for return
            Ret1.RCoef(i,:) = c;
            Ret1.PIsoMax(i) = c(1)+abs(c(2));

            % increment count to keep track of added points
            count = count +1;
            
            % Do not let program add more than 10 points
            if count >= 10 && (Ret1.PIsoMax(i) < PresMax || Ret1.BadCyc(i) == 1)


                disp(['    fit_takeuchi: Added nine points on '
                    'systolic side of curve, and Pmax']);
                disp('        remains short of actual pressure');
                disp(['        Wave ',num2str(i, '%02i'), ' is excluded']);

                Ret1.BadCyc(i) = 1;
                Ret1.VCyc(i) = 0;

                break
            end
        end

        if Ret1.VCyc(i) == true
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
if any(Ret1.VCyc)

    Ret3.ivPlotTime = [Ret3.ivPlotTime; ADD_TPoints];
    Ret3.ivPlotPres = [Ret3.ivPlotPres; ADD_PPoints];

    temp = 1:1:nfits;
    disp(['    fit_takeuchi: Vanderpool Points added on cycles ' ...
        num2str(temp(logical(Ret1.VCyc)),'%02i ')]);

    % Update ivSeg; iv1Time & iv1Pres may be updated in Vanderpool section.
    Ret2 = ivSeg;

end

% Give GUI_SINU_FIT initial ICs to work with, the average of the specific ones?
if isstruct(ICS)
    Ret1.InitIC = mean(Ret1.CycICs);
end

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(Ret1.BadCyc==1); % find indices of the bad waves
if ~isempty(indX)
    disp(['    fit_takeuchi: Some waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq(i)),'%5.3f') '.']);
    disp(['        These waves are excluded: ', num2str(indX','%02i ')]);
else
    disp(['    fit_takeuchi: All waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq(i)),'%5.3f') '.']);
end
