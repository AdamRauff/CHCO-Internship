% check to if RV_data.csv file exists
filechk = which('RV_data.csv');

% if file does not exist, creat one, and write all the headers to it
if isempty(filechk)
    fd0 = fopen('RV_data.csv', 'w');
    fprintf(fd0, 'Pnam, Pmrn, file, ');
else
    % if file exists, append to it
    fd0 = fopen('RV_data.csv', 'a');
end