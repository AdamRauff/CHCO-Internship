function [vd, vdd] = data_centdiff (dat_typ, v)
% compute higher order differences of given vector "v" for use in all 
% methods. Not totally sure if this helps, but it seems to be /slightly/
% smoother than forward diff. Also should be zero-phase lag result.

% Set timestep and data length before computing diffs
if dat_typ
    h = 1/250;
else
    h = 1/1000;
end

n = length(v);
Inv12h  = 1/(12*h);
Inv12hh = Inv12h/h;

%% Compute d(vec)/dt using 6th order central diff
% Take first three and last three terms as 4th order forward or backward
% approximations, moving evaluation point towards central.

vd(1) = Inv12h*(-25*v(1)+48*v(2)-36*v(3)+16*v(4)-3*v(5));
vd(2) = Inv12h*(-3*v(1)-10*v(2)+18*v(3)-6*v(4)+v(5));
vd(3) = Inv12h*(v(1)-8*v(2)+8*v(4)-v(5));
Inv60h = 0.2*Inv12h;
for i = 4 : 1 : n -3
    vd(i) = Inv60h*(-v(i-3)+9*v(i-2)-45*v(i-1)+45*v(i+1)-9*v(i+2)+v(i+3));
end
vd(n-2) = Inv12h*(v(n-4)-8*v(n-3)+8*v(n-1)-v(n));
vd(n-1) = Inv12h*(-v(n-4)+6*v(n-3)-18*v(n-2)+10*v(n-1)+3*v(n));
vd(n-0) = Inv12h*(3*v(n-4)-16*v(n-3)+36*v(n-2)-48*v(n-1)+25*v(n));
vd = vd';

%% Compute second derivative of vec, again 6th order cent diff
% This is for second deriviative straight from 0th order. Again here, use a
% higher order (but not quite 6th) to get more accurate edge terms.

vdd(1) = Inv12hh*(35*v(1)-104*v(2)+114*v(3)-56*v(4)+11*v(5));
vdd(2) = Inv12hh*(11*v(1)-20*v(2)+6*v(3)+4*v(4)-v(5));
vdd(3) = Inv12hh*(-v(1)+16*v(2)-30*v(3)+16*v(4)-1*v(5));
Inv180hh = Inv60h/(3*h);
for i = 4 : 1 : n -3
    vdd(i) = Inv180hh*(2*v(i-3)-27*v(i-2)+270*v(i-1)-490*v(i) ...
        +270*v(i+1)-27*v(i+2)+2*v(i+3));
end
vdd(n-2) = Inv12hh*(-v(n-4)+16*v(n-3)-30*v(n-2)+16*v(n-1)-v(n));
vdd(n-1) = Inv12hh*(-v(n-4)+4*v(n-3)+6*v(n-2)-20*v(n-1)+11*v(n));
vdd(n-0) = Inv12hh*(11*v(n-4)-56*v(n-3)+114*v(n-2)-104*v(n-1)+35*v(n));

vdd = vdd';