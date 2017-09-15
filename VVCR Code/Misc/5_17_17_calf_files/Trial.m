clear;
clc;

PathName = '';
FileName = 'csuCO214_milRV.txt';

[Pres, dPdt, Rvals, file, ~]=load_calf_p_5_17_17(PathName,FileName,100);