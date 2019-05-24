function [ivSeg, ivVal] = test_uniq_isoseg (inp_ivSeg, inp_ivVal, ivIdx, ...
    Data, isopos, runNum)

ivSeg = inp_ivSeg;
ivVal = inp_ivVal;

midl = 2;
shft = runNum - midl;

Ps = 1 + midl; Ns = Ps;
Pe = 0 + midl; Ne = Pe;

if isopos == 1
    Ps = Ps + shft;
elseif isopos == 2
    Pe = Pe - shft;
elseif isopos == 3
    Ns = Ns + shft;
elseif isopos == 4
    Ne = Ne - shft;
end

ivSeg.iv2Time(1).PosIso = ivSeg.iv2Time(1).PosIso(Ps:end);
ivSeg.iv2Pres(1).PosIso = ivSeg.iv2Pres(1).PosIso(Ps:end);
ivVal.Ps2 = ivSeg.iv2Pres(1).PosIso(1);

ivSeg.iv2Time(1).PosIso = ivSeg.iv2Time(1).PosIso(1:end-Pe);
ivSeg.iv2Pres(1).PosIso = ivSeg.iv2Pres(1).PosIso(1:end-Pe);
ivVal.Pe2 = ivSeg.iv2Pres(1).PosIso(end);

ivSeg.iv2Time(1).NegIso = ivSeg.iv2Time(1).NegIso(Ns:end);
ivSeg.iv2Pres(1).NegIso = ivSeg.iv2Pres(1).NegIso(Ns:end);
ivSeg.iv2dPdt(1).NegIso = ivSeg.iv2dPdt(1).NegIso(Ns:end);
ivVal.Ns2 = ivSeg.iv2Pres(1).NegIso(1);
ivVal.dNs2 = ivSeg.iv2dPdt(1).NegIso(1);

ivSeg.iv2Time(1).NegIso = ivSeg.iv2Time(1).NegIso(1:end-Ne);
ivSeg.iv2Pres(1).NegIso = ivSeg.iv2Pres(1).NegIso(1:end-Ne);
ivSeg.iv2dPdt(1).NegIso = ivSeg.iv2dPdt(1).NegIso(1:end-Ne);
ivVal.Ne2 = ivSeg.iv2Pres(1).NegIso(end);
ivVal.dNe2 = ivSeg.iv2dPdt(1).NegIso(end);

if isopos == 1
%   disp([ivVal.Ps2 Data.Pres_D(ivIdx.Ps1_D+runNum)])
    ivVal.dPs2 = Data.dPdt_D(ivIdx.Ps1_D+runNum);
else
%   disp([ivVal.Ps2 Data.Pres_D(ivIdx.Ps1_D+2)])
    ivVal.dPs2 = Data.dPdt_D(ivIdx.Ps1_D+2);
end

if isopos == 2
%   disp([ivVal.Pe2 Data.Pres_D(ivIdx.Pe2_D+runNum-4)])
    ivVal.dPe2 = Data.dPdt_D(ivIdx.Pe2_D+runNum-4);
else
%   disp([ivVal.Pe2 Data.Pres_D(ivIdx.Pe2_D-2)])
    ivVal.dPe2 = Data.dPdt_D(ivIdx.Pe2_D-2);
end

if isopos == 1
    fprintf ('ivSeg.Ps : ');
    for i = 1: runNum
        fprintf ('%10.6f ', 0);
    end
    for i = 1: 7-runNum
        fprintf ('%10.6f ', ivSeg.iv2Pres(1).PosIso(i));
    end
    fprintf ('\n');
elseif isopos == 2
    fprintf ('ivSeg.Pe : ');
    endst = length(ivSeg.iv2Pres(1).PosIso);
    for i = endst-runNum-2: endst
        fprintf ('%10.6f ', ivSeg.iv2Pres(1).PosIso(i));
    end
    for i = 1: 4-runNum
        fprintf ('%10.6f ', 0);
    end
    fprintf ('\n');
elseif isopos == 3
    fprintf ('ivSeg.Ns : ');
    for i = 1: runNum
        fprintf ('%10.6f ', 0);
    end
    for i = 1: 7-runNum
        fprintf ('%10.6f ', ivSeg.iv2Pres(1).NegIso(i));
    end
    fprintf ('\n');
elseif isopos == 4
    fprintf ('ivSeg.Ne : ');
    endst = length(ivSeg.iv2Pres(1).NegIso);
    for i = endst-runNum-2: endst
        fprintf ('%10.6f ', ivSeg.iv2Pres(1).NegIso(i));
    end
    for i = 1: 4-runNum
        fprintf ('%10.6f ', 0);
    end
    fprintf ('\n');
end

end

