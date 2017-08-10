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
tp = 2*pi;
P4T = P(4)*TN;

% Equation for fitting isovolumic contraction
t13 = t1-P(3);
p0m = @(P,t) a(1)/2*(P(1)-P(2))+P(2)+(P(1)-P(2))*( ...
    a(2)*cos(tp*1*t/P4T) + b(2)*sin(tp*1*t/P4T) + ...
    a(3)*cos(tp*2*t/P4T) + b(3)*sin(tp*2*t/P4T) + ...
    a(4)*cos(tp*3*t/P4T) + b(4)*sin(tp*3*t/P4T) + ...
    a(5)*cos(tp*4*t/P4T) + b(5)*sin(tp*4*t/P4T) + ...
    a(6)*cos(tp*5*t/P4T) + b(6)*sin(tp*5*t/P4T) + ...
    a(7)*cos(tp*6*t/P4T) + b(7)*sin(tp*6*t/P4T);

% First "half" of the fit residuals
zero = p0m(P,t13)-Pd;

% Time derivative of p0m, used for fitting dP/dt during isovolumic
% relaxation.
p1m = @(P,t) tp*(P(1)-P(2))/P4T*( ...
    -1*a(2)*sin(tp*1*t/P4T) + 1*b(2)*cos(tp*1*t/P4T) ...
    -2*a(3)*sin(tp*2*t/P4T) + 2*b(3)*cos(tp*2*t/P4T) ...
    -3*a(4)*sin(tp*3*t/P4T) + 3*b(4)*cos(tp*3*t/P4T) ...
    -4*a(5)*sin(tp*4*t/P4T) + 4*b(5)*cos(tp*4*t/P4T) ...
    -5*a(6)*sin(tp*5*t/P4T) + 5*b(6)*cos(tp*5*t/P4T) ...
    -6*a(7)*sin(tp*6*t/P4T) + 6*b(7)*cos(tp*6*t/P4T);

% Code to compute (dP/dt)min offset: Given the fit coefficients in P, find
% the fitted time of (dP/dt)min, then compute difference. Voila!
Tspan = [0:0.01:TN*1.2];
dPt0 = p1m (P,Tspan);
[~,idx] = min(dPt0);
tshift = Tspan(idx)-tM;

% Second "half" of the fit residuals
t23 = t2-P(3)+tshift;
zero = [zero p1m(P,t23)-dPd];

end

