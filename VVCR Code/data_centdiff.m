function [vec_ret] = data_centdiff (dat_typ, vec)
% compute higher order difference of given vector (not sure if this helps,
% but it seems to be /slightly/ smoother than forward diff. Also should be
% zero-phase lag result.

% Set timestep and data length before computing diffs
if dat_typ
    h = 1/250;
else
    h = 1/1000;
end

n = length(vec);

%% Compute d(vec)/dt using 6th order central diff
% The stuff commented out is just hunterk playing with various difference 
% To see which had the best "behavior". The 4th order forward and backward
% differences are good (low phase shift) but are still more "peaky" than
% the central diff. In the end, went with a 6th order central diff, with
% the ends approximated to best accuracy as we move in. The 4th and 6th
% order central diffs were very similar, but it is easy enough to take the
% higher order approx.

% dpf1 = diff(vec)/h;
% dpf1(n) = dpf1(n-1);
% 
% InvSixH = 1/(6*h);
% for i = 1: 1: n-3
%     dpf4(i) = InvSixH*(2*vec(i+3)-9*vec(i+2)+18*vec(i+1)-11*vec(i));
% end
% dpf4(n-2) = dpf4(n-3);
% dpf4(n-1) = dpf4(n-3);
% dpf4(n)   = dpf4(n-3);
% 
% InvTwlvH = InvSixH/2;
% for i = 3: 1: n-2
%     dpc4(i) = InvTwlvH*(-vec(i+2)+8*vec(i+1)-8*vec(i-1)+vec(i-2));
% end
% dpc4(1) = dpc4(3);
% dpc4(2) = dpc4(3);
% dpc4(n-1) = dpc4(n-2);
% dpc4(n) = dpc4(n-2);
% 
% for i = 4: 1: n
%     dpb4(i) = InvSixH*(-2*vec(i-3)+9*vec(i-2)-18*vec(i-1)+11*vec(i));
% end
% dpb4(1) = dpb4(4);
% dpb4(2) = dpb4(4);
% dpb4(3) = dpb4(4);
% 
% InvSxtyH = InvTwlvH/5;

InvSxtyH = 1/(60*h);
dpc6(1) = (vec(2)-vec(1))/h;
dpc6(2) = (vec(3)-vec(1))/(2*h);
dpc6(3) = (-vec(5)+8*vec(4)-8*vec(2)+vec(1))/(12*h);
for i = 4 : 1 : n -3
    dpc6(i) = InvSxtyH*(-vec(i-3)+9*vec(i-2)-45*vec(i-1)+45*vec(i+1)...
        -9*vec(i+2)+vec(i+3));
end
dpc6(n-2) = (-vec(n-4)+8*vec(n-3)-8*vec(n-1)+vec(n))/(12*h);
dpc6(n-1) = (vec(n-2)-vec(n))/(2*h);
dpc6(n-0) = (vec(n-1)-vec(n))/h;

% figure;
% plot(dpf1); hold on;
% plot(dpf4);
% plot(dpc4);
% plot(dpb4);
% plot(dpc6);
% pause;

%% Compute second derivative of vec, again 6th order cent diff
% This is just more experimentation; not used in "production" code. We just
% take one derivative at a time.

% ddp0 = diff(dpf1)/h;
% ddp0(n) = ddp0(n-1);

% ddpc6(1) = (dpc6(2)-dpc6(1))/h;
% ddpc6(2) = (dpc6(3)-dpc6(1))/(2*h);
% ddpc6(3) = (-dpc6(5)+8*dpc6(4)-8*dpc6(2)+dpc6(1))/(12*h);
% for i = 4 : 1 : n -3
%     ddpc6(i) = InvSxtyH*(-dpc6(i-3)+9*dpc6(i-2)-45*dpc6(i-1)+45*dpc6(i+1)...
%         -9*dpc6(i+2)+dpc6(i+3));
% end
% ddpc6(n-2) = (-dpc6(n-4)+8*dpc6(n-3)-8*dpc6(n-1)+dpc6(n))/(12*h);
% ddpc6(n-1) = (dpc6(n-2)-dpc6(n))/(2*h);
% ddpc6(n-0) = (dpc6(n-1)-dpc6(n))/h;

% figure;
% plot(ddp0); hold on;
% plot(ddpc6)
% pause;

% Set return derivative value to 6th order result.
vec_ret = dpc6';