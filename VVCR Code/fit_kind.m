function [Ret1, Ret2] = fit_kind (ivSeg, ivIdx, Data, WgtFlg)
%
% ivSeg  - Struct of all pres and time fitting values:
%            iv1Pres/iv1Time/iv2Pres/iv2Time 1st level structs; Time labels
%              actually contain indices, not real times.
%            PosIso/NegIso 2nd level structs; values for contract and relax
% Data   - Structure containing Time and Pressure.
% ICS    - Struct or Vector; If called from VVCR_ (first call), is struct
%          needed for individual-cycle ICs; if called from a GUI, contains
%          contstant initial conditions for fit.

% This causes the code to more strongly weight the isovolumic contraction
% residuals compared to the "balanced" mean weights. The value (which should be
% larger than 1) is the actual weight. If set to zero, this will have no effect.
WGHT_CONT = 0;
BOUND_VIO = 0;

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

% Ploting vectors of the fitting data for GUI_FitKind  
Ret2.iv2PlotTime = [];
Ret2.iv2PlotPres = [];
Ret2.iv2TShift = zeros(nfits,1);

% scroll through the number of rows (pressure waves) in the structures:
% ivSeg.iv2Time and ivSeg.iv2Pres
for i = 1:nfits

    % zero is computed to offset each cycle to start at t=0 (not sure why I
    % didn't do this from the start, D'OH!
    zero = Data.Time_D(ivSeg.iv2Time(i).PosIso(1));

    % Times for (dP/dt)max, (dP/dt)min, and the average period length.
    dPtimes = [Data.Time(ivIdx.dPmax2(i))-zero ...
        Data.Time(ivIdx.dPmin2(i))-zero Data.time_per];

    % Indices for the isovolumic phases. Offset of the times occurs when these
    % are used to get values from Time_D.
    posidx = ivSeg.iv2Time(i).PosIso;
    negidx = ivSeg.iv2Time(i).NegIso;

    if WgtFlg
        P0_weight = mean(ivSeg.iv2dPdt(i).NegIso)/mean(ivSeg.iv2Pres(i).PosIso);
        if WGHT_CONT
            P0_weight = P0_weight/WGHT_CONT;
        end
        nam = ' (weighted)';
    else
        P0_weight = 1;
        nam = '';
    end

%-[Examining effect of weighting
%   disp(['Cycle #' num2str(i,'%02i') ' Weighting Factor = ' ...
%       num2str(P0_weight, '%8.3f')]);

    sin_fun2 = @(P) imbedded_kind (P, Data.Time_D(posidx)-zero, ...
        Data.Time_D(negidx)-zero, dPtimes, ivSeg.iv2Pres(i).PosIso, ...
        ivSeg.iv2dPdt(i).NegIso, P0_weight);
%-[Examining time shift @IC and at finish
%       ivSeg.iv2dPdt(i).NegIso, 0);
    
    % Deriving the initial values from the data
    % P1 2.5*Pes - doesn't depend on any other fit, and is a good start (this
    %    value gives a VVCR of 1.5, halfway between conditions of max work (1.0)
    %    and max efficiency (2.0).
    % P2 Pmin, guess small, like 10.
    % P3 use 0.0 (start of time-normalized iso contraction)
    % P4 use 58% (that's their IC!)

    j = ivIdx.goodcyc2(i);
    c2 = [2.5*Data.PesP(j), 10, 0.00, 0.58];
    Ret1.CycICs(i,:) = c2; 

    % First Set of Limits - very weak bounds on t_Pmax. OLD DO NOT USE
    % lb = [Data.Pes2(i)  0.0            -0.1   0.2];
    % ub = [        1000 30.0  Data.Time(end)   0.8];
    %
    % New Limits:
    % - Max pressure must be larger than Pes, but have a "reasonable" upper
    %   bound (500 mmHg?)
    % - Min pressure must be larger than zero, but less than... Pes? Or
    %   something even smaller? 
    % - t0 should be pretty small and probably positive
    % - Beta: 0.6 is the "best value" for rats, so we give it some leeway.

    lb = [Data.PesP(i)  0.0 -0.005 0.48];
    ub = [         500 40.0  0.020 0.72];

    [c,SSE,~] = lsqnonlin (sin_fun2,c2,lb,ub,opts1);
    
    % r^2 value; if the fit was bad, mark that wave.
    WavePs = [ivSeg.iv2Pres(i).PosIso; ivSeg.iv2dPdt(i).NegIso];
    SSTO = norm(WavePs-mean(WavePs))^2;
    Ret1.Rsq(i) = 1-SSE/SSTO;
    
%-[Examining effect of weighting
%   disp(['    Rsq ' num2str(Ret1.Rsq(i), '%6.4f')]);
       
    if Ret1.Rsq(i) < 0.60
       Ret1.BadCyc(i) = 1; 
    end

    if any( abs(c-lb) < 1e-6 ) || any ( abs(ub-c) < 1e-6 )
        if ~BOUND_VIO 
            fprintf('    fit_kind%s: fit bounds violated on cycle %02i', ...
                nam, i);
            BOUND_VIO = 1;
        else
            fprintf(' %02i', i);
        end    
        Ret1.BadCyc(i) = 1;
    end
    
    %getting all the c values in a matrix
    Ret1.RCoef(i,:) = c; 

%   fprintf ('        RevPmax %7.3f %6.3f %6.4f %6.4f\n', c(1), c(2), c(3), c(4));
    
    % store the time points and pressure points in one array for easy plotting -
    % first pass (call from VVCR_); otherwise, reconsitute these arrays if
    % needed just outside this loop.
    [~, tsh, padd] = data_kind (c, Data.Time_D(posidx(1))-zero, dPtimes);
    psh = padd-Data.Pres_D(ivIdx.dPmin2_D(i));

    Ret2.iv2TShift(i) = tsh;
    
    Ret2.iv2PlotTime = [Ret2.iv2PlotTime Data.Time_D(posidx) ...
        tsh+Data.Time_D(negidx)];
    Ret2.iv2PlotPres = [Ret2.iv2PlotPres Data.Pres_D(posidx)' ...
        psh+Data.Pres_D(negidx)'];

%-[Examining time shift @IC and at finish
% [~, ICshift] = imbedded_kind (c2, Data.Time_D(posidx), ...
%     Data.Time_D(negidx), dPtimes, 0, 0, 1);
% [~, FTshift] = imbedded_kind (c, Data.Time_D(posidx), ...
%     Data.Time_D(negidx), dPtimes, 0, 0, 1);
% disp(['Fit # ' num2str(i, '%02i')]);
% disp(['    ICs Shift ', num2str(ICshift, '%6.4f ') ' diff ' ...
%     num2str(ICshift(1)-ICshift(2), ' %+6.4f') ])
% disp(['    Fit Shift ', num2str(FTshift, '%6.4f ') ' diff ' ...
%     num2str(FTshift(1)-FTshift(2), ' %+6.4f') ])
% disp(['    Max Pressure ', num2str(c(1), '%8.3f') ' Ps1 = ' ...
%     num2str(Data.Time_D(ivIdx.Ps1_D(i)), '%6.4f') ' C(3,4) = ' ...
%     num2str(c(3:4), '%6.4f ') ' ' num2str(Data.Time_D(ivIdx.Ps1_D(i))- ...
%     c(3),'(%+6.4f)') ]);

end

% print to command line the waves that were not fit correctly. This is used as a
% debugger to check that the "bad" waves, the ones that don't have a good fit,
% are not utilized in the VVCR calculation.
indX = find(Ret1.BadCyc==1); % find indices of the bad waves
if BOUND_VIO
    fprintf('\n');
end
if ~isempty(indX)
    disp(['    fit_kind' nam ': Some waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq(Ret1.BadCyc~=1)),'%5.3f') '.']);
    disp(['        These waves are excluded: ', num2str(indX','%02i ')]);
else
    disp(['    fit_kind' nam ': All waves fit well, ave R^2 = ' ...
            num2str(mean(Ret1.Rsq),'%5.3f') '.']);
end

%fprintf ('        RevPmax %7.3f %6.3f %6.4f %6.4f\n', c(1), c(2), c(3), c(4));

% END OF fit_kind
end

function [ zero ] = imbedded_kind ( P, t1, t2, tM, Pd, dPd, weight )
%IMBEDDED_KIND Summary of this function goes here
%   This is a placeholder function with multiple arguments that will be
%   passed using an anoymous handle to lsqnonlin fit within VVCR_FINAL_*
%   (more explanation to come in that script). lsqnonlin requires an
%   function with a single argument (namely, the fit coefficients P)
%     Input Arguments:
%       P      - fit coefficients Pmax, Pmin, t0, tpmax
%       t1     - time vector for isovolumic contraction (fit to P)
%       t2     - time vector for isovolumic relaxation (fit to dP/dt)
%       tM     - times of (dP/dt)max, (dP/dt)min in actual data, and actual
%                period length
%       Pd     - pressure data during isovolumic contraction
%       dPd    - dP/dt data during isovolumic relaxation
%       weight - Error weighting for Pd error (= mean(Pd)/mean(dPd))
%     Output Arguments:
%       zero - difference between fit and data (combined vector of 
%              p0m-Pd and p1m-dPd)
%
%-[Examining time shift @IC and at finish
%%function [ varargout ] = imbedded_kind ( P, t1, t2, tM, Pd, dPd, flag )

% Coefficients for the multiharmonic fit
a = [1.0481 -0.4361 -0.0804 -0.0148  0.0020  0.0023  0.0012];
b = [ 0.000  0.2420 -0.0255 -0.0286 -0.0121 -0.0039 -0.0016];
TN = 2.658;

% Substitutions to simplify eqns and reduce math overhead. t_pmax is the
% (normalized) time of (dP/dt)max plus beta*(diff between extrema); t0 [P(3)] is
% not subtracted from this because each cycle time is iso-normalized such that
% the first point of the isovolumic contraction occurs at t_N = 0. Thus the
% normalized position of the isovolumic max is only relative to beta [P(4)].
% Then tp_P2T is the constant in the sin/cos terms.
t_pmax = tM(1) + P(4)*(tM(2)-tM(1));
tp_P2T = 2*pi/(t_pmax*TN);

% Equation for fitting isovolumic contraction: straight out of the paper. t13 is
% the (iso-normalized) time minus t0 (fitting constant P(3)).
t13 = t1'-P(3);
p0m = @(P,t) a(1)/2*(P(1)-P(2))+P(2)+(P(1)-P(2))*( ...
    a(2)*cos(tp_P2T*1*t) + b(2)*sin(tp_P2T*1*t) + ...
    a(3)*cos(tp_P2T*2*t) + b(3)*sin(tp_P2T*2*t) + ...
    a(4)*cos(tp_P2T*3*t) + b(4)*sin(tp_P2T*3*t) + ...
    a(5)*cos(tp_P2T*4*t) + b(5)*sin(tp_P2T*4*t) + ...
    a(6)*cos(tp_P2T*5*t) + b(6)*sin(tp_P2T*5*t) + ...
    a(7)*cos(tp_P2T*6*t) + b(7)*sin(tp_P2T*6*t) );

% First "half" of the fit residuals.
zero = weight*(p0m(P,t13)-Pd);

% Time derivative of p0m, used for fitting dP/dt during isovolumic relaxation.
p1m = @(P,t) tp_P2T*(P(1)-P(2))*( ...
    -1*a(2)*sin(tp_P2T*1*t) + 1*b(2)*cos(tp_P2T*1*t) ...
    -2*a(3)*sin(tp_P2T*2*t) + 2*b(3)*cos(tp_P2T*2*t) ...
    -3*a(4)*sin(tp_P2T*3*t) + 3*b(4)*cos(tp_P2T*3*t) ...
    -4*a(5)*sin(tp_P2T*4*t) + 4*b(5)*cos(tp_P2T*4*t) ...
    -5*a(6)*sin(tp_P2T*5*t) + 5*b(6)*cos(tp_P2T*5*t) ...
    -6*a(7)*sin(tp_P2T*6*t) + 6*b(7)*cos(tp_P2T*6*t) );

% Code to compute (dP/dt)min offset: Given the fit coefficients in P, find the
% fitted time of (dP/dt)min, then compute difference. Voila!
%
% In more detail: P coefficients determine point of maximum for multiharmonic
% fit. So we just compute that. Then, we also already know the time at which the
% actual (dP/dt)min occurs. These are independent events, given a specific set
% of P coefficients. Then, knowing both of these times, we can choose the time
% for the multiharmonic at which comparisons are made - and it's centered around
% each vector's (dP/dt)min.
%
% Tspan is just some range (shifted or not it's irrelevant) in which we believe
% we can find the isovolumic (dP/dt)min.
Tspan = t1(end) : 0.0001: tM(3)*0.7;
dPt0 = p1m (P, Tspan);
[~,idx] = min(dPt0);

% Once isovolumic (dP/dt)min is found in [Tspan(idx)] we find the shift distance
% to the minimum found in the data [tM(2)]. Because this match will be negative
% shifted by t0 [P(3)] in the function, we add it back here.
tshift = Tspan(idx)-(tM(2)-P(3));

%-[Examining time shift @IC and at finish
%if flag
%    varargout{2} = [Tspan(idx) tM(2)-P(3)];
%    varargout(3:nargout) = {[]};
%else
%    varargout(2:nargout) = {[]};
%end

% Second "half" of the fit residuals
t23 = t2'-P(3)+tshift;

%-[Residuals in each section]
%maxr1 = norm(zero);
%maxr2 = norm(p1m(P,t23)-dPd);
%disp(['    Norm Residuals ' num2str(maxr1) ' ' num2str(maxr2)]);

zero = [zero; (p1m(P,t23)-dPd)];

end
