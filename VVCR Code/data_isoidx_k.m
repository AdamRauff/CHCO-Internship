function [ivIdx, ivVal, badcyc] = data_isoidx (mysz, Dat, Ext, ivIdx, ...
    ivVal, badcyc)
% Find isovolumic timings for Kind method points.
%
% Kind method states that th pos iso starts at R wave, which isn't available,
% so we approximate by using the same time as Takeuchi. However, when waves
% are deleted we need unique vectors for each method, so we copy the Takeuchi
% (iv*.Ps1) vetors into the Kind vectors. We also create unique copies of the
% extrema vectors for the same reason.

disp('    data_isoidx: finding Kind indices');

ivIdx.Ps2 = ivIdx.Ps1;
ivVal.Ps2 = ivVal.Ps1; 
ivVal.Pe2 = zeros(mysz,1);
ivIdx.Pe2 = zeros(mysz,1);
ivVal.Ns2 = zeros(mysz,1);
ivIdx.Ns2 = zeros(mysz,1);
ivVal.Ne2 = zeros(mysz,1);
ivIdx.Ne2 = zeros(mysz,1);

ivIdx.dPmax2 = Ext.dPmaxIdx;
ivVal.dPmax2 = Ext.dPmaxVal;
ivIdx.dPmin2 = Ext.dPminIdx;
ivVal.dPmin2 = Ext.dPminVal;

% I guess originally I just avoided using the first cycle if it was bad,
% but really I must avoid using any Ps1 values that are bad - those are the
% positive badcyc.T values.
[~,idx] = find(badcyc.T > 0);
if isempty(idx)
    badcyc.K = [];
else
%   badcyc.K = 1;
    badcyc.K = badcyc.T(idx);
end

% # of points it can look before it stops... each point is 4ms.
lenlim = 10;

for i = 1:mysz
    % Find pres at (dP/dt)max, then find 10% increase, this is end of Kind
    % pos iso segment. 
    Pcut = 1.10*Dat.Pres(ivIdx.dPmax2(i));
    Pend = ivIdx.dPmax2(i);
    while Dat.Pres(Pend) < Pcut
        Pend = Pend + 1;

        % Don't let it leave this cycle. If we don't see a 10% rise from
        % (dP/dt)max, simply take the maximum pressure as the end.
        if Pend > ivIdx.dPmin2(i)
            [~, Pend] = max(Dat.Pres(ivIdx.dPmax2(i):ivIdx.dPmin2(i)));
            break;
        end
    end
    ivVal.Pe2(i) = Dat.Pres(Pend);
    ivIdx.Pe2(i) = Pend;

    % Find pres at (dP/dt)min, then find ±20% change, these are start and end
    % of Kind neg iso segment.
    Pcut = 1.20*Dat.Pres(ivIdx.dPmin2(i));
    Pstr = ivIdx.dPmin2(i);
    Plim = Pstr - lenlim;
    while Dat.Pres(Pstr) < Pcut
        Pstr = Pstr - 1;

        % Don't let this go too far: add this cycle to "bad list" if we back up
        % to (dP/dt)max (which is really far!) or if we go farther than 5 points
        % away (this may need revision).
        if Pstr < ivIdx.dPmax2(i) || Pstr == Plim
            disp(['        curve # ' num2str(i, '%02i') ', can''t find ' ...
                'start of negiso segment, skipping.']);
            badcyc.K = [badcyc.K, i];
            break;
        end

    end
    ivVal.Ns2(i) = Dat.Pres(Pstr);
    ivIdx.Ns2(i) = Pstr;

    Pcut = 0.80*Dat.Pres(ivIdx.dPmin2(i));
    Pend = ivIdx.dPmin2(i);
    Dend = length(Dat.Pres);
    Plim = Pend + lenlim;
    while Dat.Pres(Pend) > Pcut
        Pend = Pend + 1;

        % Don't let this go too far: add this cycle to "bad list" if we exceed
        % the Takeuchi end (i.e. go beyond the isovolumic portion), we go
        % farther than 5 points away (this may need revision), or, if we're on
        % the last cycle, don't go off the end of the data.
        if Pend > Dend
            disp(['        curve # ' num2str(i, '%02i') ', end of negiso ' ...
	        'segment not captured at end of sample, skipping.']);

            Pend = Dend;
            badcyc.K = [badcyc.K, i];
            break;
        elseif Pend >= ivIdx.Ne1(i) || Pend == Plim
            disp(['        curve # ' num2str(i, '%02i') ', can''t find '...
                'end of negiso segment, skipping.']);

            badcyc.K = [badcyc.K, i];
            break;
        end

    end
    ivVal.Ne2(i) = Dat.Pres(Pend);
    ivIdx.Ne2(i) = Pend;

end
