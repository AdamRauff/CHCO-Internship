function [ zero ] = pmax_multiharm( P, t1, t2, tM, TN, Pd, dPd )
%PMAX_MULTIHARM Summary of this function goes here
%   This is a placeholder function with multiple arguments that will be
%   passed using an anoymous handle to lsqnonlin fit within VVCR_FINAL_*
%   (more explanation to come in that script). lsqnonlin requires an
%   function with a single argument (namely, the fit coefficients P)
%       Input Arguments:
%           P   - fit coefficients Pmax, Pmin, t0, tpmax
%           t1  - time vector for isovolumic contraction (fit to P)
%           t2  - time vector for isovolumic relaxation (fit to dP/dt)
%           tM  - time of dP/dt in actual data
%           TN  - period of actual data
%           Pd  - pressure data during isovolumic contraction
%           dPd - dP/dt data during isovolumic relaxation
%       Output Arguments:
%           zero - difference between fit and data (combined vector of 
%                  p0m-Pd and p1m-dPd)
%

% Coefficients for the multiharmonic fit
a = [10.481 -0.4361 -0.0804 -0.0148  0.0020  0.0023  0.0012]
b = [ 0.000  0.2420 -0.0255 -0.0286 -0.0121 -0.0039 -0.0016]

% Substitutions to simplify eqns and reduce math overhead
tp_P2T = 2*pi/(P(4)*TN);

% Equation for fitting isovolumic contraction: straight out of the paper.
t13 = t1-P(3);
p0m = @(P,t) a(1)/2*(P(1)-P(2))+P(2)+(P(1)-P(2))*( ...
    a(2)*cos(tp_P2T*1*t) + b(2)*sin(tp_P2T*1*t) + ...
    a(3)*cos(tp_P2T*2*t) + b(3)*sin(tp_P2T*2*t) + ...
    a(4)*cos(tp_P2T*3*t) + b(4)*sin(tp_P2T*3*t) + ...
    a(5)*cos(tp_P2T*4*t) + b(5)*sin(tp_P2T*4*t) + ...
    a(6)*cos(tp_P2T*5*t) + b(6)*sin(tp_P2T*5*t) + ...
    a(7)*cos(tp_P2T*6*t) + b(7)*sin(tp_P2T*6*t);

% First "half" of the fit residuals
zero = p0m(P,t13)-Pd;

% Time derivative of p0m, used for fitting dP/dt during isovolumic relaxation.
p1m = @(P,t) tp_P2T*(P(1)-P(2))*( ...
    -1*a(2)*sin(tp_P2T*1*t) + 1*b(2)*cos(tp_P2T*1*t) ...
    -2*a(3)*sin(tp_P2T*2*t) + 2*b(3)*cos(tp_P2T*2*t) ...
    -3*a(4)*sin(tp_P2T*3*t) + 3*b(4)*cos(tp_P2T*3*t) ...
    -4*a(5)*sin(tp_P2T*4*t) + 4*b(5)*cos(tp_P2T*4*t) ...
    -5*a(6)*sin(tp_P2T*5*t) + 5*b(6)*cos(tp_P2T*5*t) ...
    -6*a(7)*sin(tp_P2T*6*t) + 6*b(7)*cos(tp_P2T*6*t);

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
Tspan = [0:0.01:TN*1.2];
dPt0 = p1m (P,Tspan);
[~,idx] = min(dPt0);
tshift = Tspan(idx)-tM;

% Second "half" of the fit residuals
t23 = t2-P(3)+tshift;
zero = [zero p1m(P,t23)-dPd];

end
