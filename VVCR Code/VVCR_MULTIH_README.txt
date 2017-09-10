(1) reads data from appropriate file type based on name or extension then check for data validity
(2) Filter data (if required) and construct time array
(3) Using findpeaks, determine dP/dt minima, maxima on whole set; check these max/min sets for completeness (no min before first max, no max after final min)
(4) Call the "cycle check" GUI for check of max/min gating prior to running fit (FIRST GUI APPEARS)
    GUI return processing similar to after 2nd GUI
(5) Adjusted peaks back from GUI are used to find EDP (start of isovolumic contraction) in turn end of isovolumic relaxation
(6) Remove bad curves after finding these isovolumic points
(7) Create the 2x oversampled dataset
(8) Fit the data
(9) Call the "fit check" GUI to show fitted isovolumic curves
(10) pack data as appropriate for return to runAll
    GUI return processing similar to after 1st GUI

Variable Name Mappings are in subst (perl script)

Created in data_filter:
Data.time_step
Data.time_end
timeOld/Data.Time
[time/Data.Time]
Data.dPdt
Data.Pres
Data.dPdtAlt
Data.PresAlt

Created in data_maxmin: (originally Extrema but later, ivIdx/ivVal)
dPmaxidx/ivIdx.dPmax
ivVal.dPmax
dPminidx/ivIdx.dPmin
ivVal.dPmin

Created in data_isoidx:
iv1PsIdx/ivIdx.Ps1
iv1PsVal/ivVal.Ps1
iv1NeIdx/ivIdx.Ne1
ivVal.Ne1
ivIdx.Pe2
ivVal.Pe2
ivIdx.Ns2
ivVal.Ns2
ivIdx.Ne2
ivVal.Ne2

Created in data_double:
timeDoub/Data.Time_D
PresDoub/Data.Pres_D
Data.P_es
Data.P_esTimes

Created in data_isoseg:
iv1PsIdx_Doub / .i1d -> (ivIdx).Ps1_D
iv1NeIdx_Doub / .i2d -> (ivIdx).Ne1_D
ivIdx.Pe2_D
ivIdx.Ns2_D
ivIdx.Ne2_D

isoFitPts/ivSeg
ivSeg.iv1Time.PosIso
ivSeg.iv1Time.NegIso
ivSeg.iv1Pres.PosIso
ivSeg.iv1Pres.NegIso
ivSeg.iv2Time.PosIso
ivSeg.iv2Time.NegIso
ivSeg.iv2Pres.PosIso
ivSeg.iv2Pres.NegIso


