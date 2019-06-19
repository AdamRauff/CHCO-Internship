function [ivVal, ivIdx] = clean_isoidx (badcyc, ivVal, ivIdx, meth, field)
% Generalized cleaning of bad curves out of the value and index vectors for
% VVCR landmarks using field names passed in the string variable 'fields'.
%
% Relevant field names are: ['PsX  '; 'PeX  '; 'NsX  '; 'NeX  '], where
% X is a integer from 1-3, ['dPmaxX' 'dPminX'] where again X is the same 
% integer, and 'Pes3'.

nam = sprintf ('%s%1i', 'goodcyc', meth);
ivIdx.(nam) = setdiff(1:1:ivIdx.mxIdx, badcyc);

kk = size(field);
kk = kk(1);

if ~isempty(badcyc)
    for i = length(badcyc):-1:1
        j = badcyc(i);
        for k = 1:1:kk  
            nam = strtrim(field(k,:));
            ivVal.(nam)(j) = [];
            ivIdx.(nam)(j) = [];
        end
    end
end
