function [Ret1, Ret2] = fit_kind (ivSeg, ivIdx, Data, MeanTPmax, WgtFlg)
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
% residuals compared to the "balanced" mean weights. The value (which
% should be larger than 1) is the actual weight. If set to zero, this will
% have no effect.
WGHT_CONT = 0;

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

% scroll through the number of rows (pressure waves) in the
% structures: ivSeg.iv2Time and ivSeg.iv2Pres
for i = 1:nfits

    % Times for (dP/dt)max, (dP/dt)min, and the average period length
    dPtimes = [Data.Time(ivIdx.dPmax2(i)) Data.Time(ivIdx.dPmin2(i)) ...
        Data.time_per];

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

    sin_fun2 = @(P) imbedded_kind (P, Data.Time_D(posidx), ...
        Data.Time_D(negidx), dPtimes, ivSeg.iv2Pres(i).PosIso, ...
        ivSeg.iv2dPdt(i).NegIso, P0_weight);
%-[Examining time shift @IC and at finish
%       ivSeg.iv2dPdt(i).NegIso, 0);

    % Start of isovolumic contraction, isovolumic contration duration. Used
    % in ICs and fitting limits.
    tstivc = Data.Time_D(ivIdx.Ps2_D(i));
    ivcdur = Data.Time_D(ivIdx.Pe2_D(i))-tstivc;
    
    % Deriving the initial values from the data
    % P1 Mean Pmax from Takeuchi fits.
    % P2 Pmin, guess small, like 2-5.
    % P3 use Time(ivIdx.Ps1) - that will be close
    % P4 use 58% (that's their IC!)

    c2 = [MeanTPmax, 2, tstivc, 0.58];
    Ret1.CycICs(i,:) = c2; 

    % First Set of Limits - very weak bounds on t_Pmax.
    % lb = [Data.Pes2(i)  0.0            -0.1   0.2];
    % ub = [        1000 30.0  Data.Time(end)   0.8];
    %
    % New Limits:
    % - Max pressure must be larger than Pes, but have a "reasonable" upper
    %   bound (500 mmHg?)
    % - Min pressure must be larger than zero, but less than... Pes? Or
    %   something even smaller? 
    % - t_Pmax... this is the hard one. Should the fit be allowed to push
    %   the start time beyond the start of the isovolumic phase? I think no.
    % - Beta: 0.6 is the "best value" for rats, so we give it some leeway.

    lb = [Data.Pes2(i)  0.0 tstivc-ivcdur 0.50];
    ub = [         500 30.0 tstivc        0.70];

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

    if any( (c-lb) < 0 ) || any ( (ub-c) < 0 )
       disp(['    fit_kind' nam ': fit bounds violated on cycle ' ...
           num2str(i, '%02i')]);
       Ret1.BadCyc(i) = 1;
    end
    
    %getting all the c values in a matrix
    Ret1.RCoef(i,:) = c; 

    % store the time points and pressure points in one array for easy
    % plotting - first pass (call from VVCR_); otherwise, reconsitute these
    % arrays if needed just outside this loop.

    [~, tsh, padd] = data_kind (c, Data.Time_D(posidx(1)), dPtimes);
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

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(Ret1.BadCyc==1); % find indices of the bad waves
if ~isempty(indX)
    disp(['    fit_kind' nam ': Some waves fit well, ave R^2 = ' ...
        num2str(mean(Ret1.Rsq(Ret1.BadCyc~=1)),'%5.3f') '.']);
    disp(['        These waves are excluded: ', num2str(indX','%02i ')]);
else
    disp(['    fit_kind' nam ': All waves fit well, ave R^2 = ' ...
            num2str(mean(Ret1.Rsq),'%5.3f') '.']);
end

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
zero = weight*(p0m(P,t13)-Pd);

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
