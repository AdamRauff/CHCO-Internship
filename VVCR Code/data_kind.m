function [ pres, tshift, pshift ] = data_kind ( P, tdat, tM )
%   This is an external function that provides the Kind fit for plotting
%   and other non-fitting uses.
%     Input Arguments:
%       P    - fit coefficients Pmax, Pmin, t0, tpmax
%       tdat - time vector for pressure output
%       tM   - times of (dP/dt)max, (dP/dt)min in actual data, and actual
%              period length
%     Output Arguments:
%       pres   - Pressure for plotting in GUI_FitKind for check of fit.
%       tshift - Time shift to align actual data points in isovolumic
%                contraction phase with fit pressure (shift abscissa).
%       pshift - Actual Kind pressure @tshift (again, to be used to shift
%                the ordinate of the data).
%

% Coefficients for the multiharmonic fit
a = [1.0481 -0.4361 -0.0804 -0.0148  0.0020  0.0023  0.0012];
b = [ 0.000  0.2420 -0.0255 -0.0286 -0.0121 -0.0039 -0.0016];
TN = 2.658;

% Substitutions to simplify eqns and reduce math overhead
t_pmax = tM(1) + P(4)*(tM(2)-tM(1)) - P(3);
tp_P2T = 2*pi/(t_pmax*TN);

% Equation for fitting isovolumic contraction: straight out of the paper.
tdatshift = tdat'-P(3);
p0m = @(P,t) a(1)/2*(P(1)-P(2))+P(2)+(P(1)-P(2))*( ...
    a(2)*cos(tp_P2T*1*t) + b(2)*sin(tp_P2T*1*t) + ...
    a(3)*cos(tp_P2T*2*t) + b(3)*sin(tp_P2T*2*t) + ...
    a(4)*cos(tp_P2T*3*t) + b(4)*sin(tp_P2T*3*t) + ...
    a(5)*cos(tp_P2T*4*t) + b(5)*sin(tp_P2T*4*t) + ...
    a(6)*cos(tp_P2T*5*t) + b(6)*sin(tp_P2T*5*t) + ...
    a(7)*cos(tp_P2T*6*t) + b(7)*sin(tp_P2T*6*t) );

% Pressure for plotting against fitting data
pres = p0m (P, tdatshift);

if length(tdat) > 1
    tshift = 0;
    pshift = 0;
    return;
end

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
Tspan = tdatshift(1) : 0.005 : (tdatshift(1)+tM(3)*1.1);
dPt0 = p1m (P, Tspan);
[~,idx] = min(dPt0);

% These are offsets that will be used to map the isovolumic relaxation data
% points onto the Kind curve (hopefully?!?).
tshift = Tspan(idx)-(tM(2)-P(3));
pshift = p0m (P, Tspan(idx));

end
