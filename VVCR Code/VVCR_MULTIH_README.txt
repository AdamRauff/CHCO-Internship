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

