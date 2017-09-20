function [Ret1] = fit_kind (ivSeg, ivIdx, Data, FitT)
%
% ivSeg  - Struct of all pres and time fitting values:
%            iv1Pres/iv1Time/iv2Pres/iv2Time 1st level structs; Time labels
%              actually contain indices, not real times.
%            PosIso/NegIso 2nd level structs; values for contract and relax
% Data   - Structure containing Time and Pressure.
% ICS    - Struct or Vector; If called from VVCR_ (first call), is struct
%          needed for individual-cycle ICs; if called from a GUI, contains
%          contstant initial conditions for fit.

opts1 = optimoptions (@lsqnonlin);
opts1.Display = 'off';
opts1.MaxFunctionEvaluations = 2000;
opts1.MaxIterations = 1000;
opts1.FiniteDifferenceType = 'central';

% Variables for main fit, all returned in Ret1
nfits = length(ivSeg.iv2Pres);

Ret1.Rsq = zeros(nfits,1);    % Goodness of fit coefficients
Ret1.RCoef = zeros(nfits,4);  % Fit regression constants
Ret1.CycICs = zeros(nfits,4); % Which waveforms had a bad fit
Ret1.BadCyc = zeros(nfits,1); % Saved cycle-specific ICS

% scroll through the number of rows (pressure waves) in the
% structures: ivSeg.iv2Time and ivSeg.iv2Pres
for i = 1:nfits

    % Times for (dP/dt)max, (dP/dt)min, and the average period length
    dPtimes = [Data.Time(ivIdx.dPmax(i)) Data.Time(ivIdx.dPmin(i)) ...
        Data.time_per];

    sin_fun2 = @(P) pmax_multiharm (P, Data.Time_D(ivSeg.iv2Time(i).PosIso),...
        Data.Time_D(ivSeg.iv2Time(i).NegIso), dPtimes, ...
	ivSeg.iv2Pres(i).PosIso, ivSeg.iv2Pres(i).NegIso);	
    
    % Deriving the initial values from the data
    % P1 Pmax from Takeuchi(?)
    % P2 Pmin, guess small, like 2-5.
    % P3 use Time(ivIdx.Ps1) - that will be close
    % P4 use 60% (that's what's in their figure!)
    c2 = [FitT.PIsoMax(i), 2, Data.Time_D(ivIdx.Ps1_D(i)), 0.6];
    Ret1.CycICs(i,:) = c2; 

    lb = [Data.P_es(i)  0.0            -0.1   0.2];
    ub = [        1000 30.0  Data.Time(end)   0.8];
    [c,SSE,~] = lsqnonlin (sin_fun2,c2,lb,ub,opts1);
    
    % r^2 value; if the fit was bad, mark that wave.
    WavePs = [ivSeg.iv2Pres(i).PosIso; ivSeg.iv2Pres(i).NegIso];
    SSTO = norm(WavePs-mean(WavePs))^2;
    Ret1.Rsq(i) = 1-SSE/SSTO;
    
    if Ret1.Rsq(i) < 0.90
       Ret1.BadCyc(i) = 1; 
    end

    if any( (c-lb) < 0 ) || any ( (ub-c) < 0 )
       disp(['    fit_kind: fit bounds violated on cycle ' ...
           num2str(i, '%02i')]);
       Ret1.BadCyc(i) = 1;
    end

    
    %getting all the c values in a matrix
    Ret1.RCoef(i,:) = c; 
    
end

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(Ret1.BadCyc==1); % find indices of the bad waves
if ~isempty(indX)
    disp(['    fit_kind: Some waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq(Ret1.BadCyc~=1)),'%5.3f') '.']);
    disp(['        These waves are excluded: ', num2str(indX','%02i ')]);
else
    disp(['    fit_kind: All waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq),'%5.3f') '.']);
end

% END OF fit_kind
end

function [ zero ] = pmax_multiharm( P, t1, t2, tM, Pd, dPd )
%PMAX_MULTIHARM Summary of this function goes here
%   This is a placeholder function with multiple arguments that will be
%   passed using an anoymous handle to lsqnonlin fit within VVCR_FINAL_*
%   (more explanation to come in that script). lsqnonlin requires an
%   function with a single argument (namely, the fit coefficients P)
%       Input Arguments:
%           P   - fit coefficients Pmax, Pmin, t0, tpmax
%           t1  - time vector for isovolumic contraction (fit to P)
%           t2  - time vector for isovolumic relaxation (fit to dP/dt)
%           tM  - times of (dP/dt)max, (dP/dt)min in actual data, and actual
%                 period length
%           Pd  - pressure data during isovolumic contraction
%           dPd - dP/dt data during isovolumic relaxation
%       Output Arguments:
%           zero - difference between fit and data (combined vector of 
%                  p0m-Pd and p1m-dPd)
%

% Coefficients for the multiharmonic fit
a = [1.0481 -0.4361 -0.0804 -0.0148  0.0020  0.0023  0.0012];
b = [ 0.000  0.2420 -0.0255 -0.0286 -0.0121 -0.0039 -0.0016];
TN = 2.658;

% Substitutions to simplify eqns and reduce math overhead
t_pmax = tM(1) + P(4)*(tM(2)-tM(1)) - P(3);
tp_P2T = 2*pi/(t_pmax*TN);

% Equation for fitting isovolumic contraction: straight out of the paper.
t13 = t1'-P(3);
p0m = @(P,t) a(1)/2*(P(1)-P(2))+P(2)+(P(1)-P(2))*( ...
    a(2)*cos(tp_P2T*1*t) + b(2)*sin(tp_P2T*1*t) + ...
    a(3)*cos(tp_P2T*2*t) + b(3)*sin(tp_P2T*2*t) + ...
    a(4)*cos(tp_P2T*3*t) + b(4)*sin(tp_P2T*3*t) + ...
    a(5)*cos(tp_P2T*4*t) + b(5)*sin(tp_P2T*4*t) + ...
    a(6)*cos(tp_P2T*5*t) + b(6)*sin(tp_P2T*5*t) + ...
    a(7)*cos(tp_P2T*6*t) + b(7)*sin(tp_P2T*6*t) );

% First "half" of the fit residuals
zero = (p0m(P,t13)-Pd);

% Time derivative of p0m, used for fitting dP/dt during isovolumic relaxation.
p1m = @(P,t) tp_P2T*(P(1)-P(2))*( ...
    -1*a(2)*sin(tp_P2T*1*t) + 1*b(2)*cos(tp_P2T*1*t) ...
    -2*a(3)*sin(tp_P2T*2*t) + 2*b(3)*cos(tp_P2T*2*t) ...
    -3*a(4)*sin(tp_P2T*3*t) + 3*b(4)*cos(tp_P2T*3*t) ...
    -4*a(5)*sin(tp_P2T*4*t) + 4*b(5)*cos(tp_P2T*4*t) ...
    -5*a(6)*sin(tp_P2T*5*t) + 5*b(6)*cos(tp_P2T*5*t) ...
    -6*a(7)*sin(tp_P2T*6*t) + 6*b(7)*cos(tp_P2T*6*t) );

% Code to compute (dP/dt)min offset: Given the fit coefficients in P, find
% the fitted time of (dP/dt)min, then compute difference. Voila!
%
% In more detail: P coefficients determine point of maximum for multiharmonic
% fit. So we just compute that. Then, we also already know the time at which
% the actual (dP/dt)min occurs. These are independent events, given a specific
% set of P coefficients. Then, knowing both of these times, we can choose the
% time for the multiharmonic at which comparisons are made - and it's centered 
% around each vector's (dP/dt)min.
%
Tspan = t13(1) : 0.005: (t13(1)+tM(3)*1.1);
dPt0 = p1m (P, Tspan);
[~,idx] = min(dPt0);
tshift = Tspan(idx)-(tM(2)-P(3));

% Second "half" of the fit residuals
t23 = t2'-P(3)+tshift;
zero = [zero; (p1m(P,t23)-dPd)];

end
