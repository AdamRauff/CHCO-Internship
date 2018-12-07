function [Pres, dPdt, Rvals, hdat] = loadp_10_04_18(npath, file)

debug = 0;

%% Open the passed filename.
filename = strcat(npath, file);   % Get the full filename
redname = basename(filename, 3);  % Get its basename for reporting
fd0 = fopen(filename);            % open the file descriptor
if fd0 < 1
    [Pres, dPdt, Rvals, hdat] = deal (-1);
    return;
end

%% Read in entire file to cell
all = textscan(fd0, '%s', 'Delimiter', '\n');
all = all{1,1};

retval = fclose(fd0);
if retval ~= 0
    [Pres, dPdt, Rvals, hdat] = deal (-1);
    return;
end

%% Process Whatever Exists in Header
header = 1;
while isempty(regexp(all{header},'Begin'))
    header = header+1;
end

hdat = struct('Nam', '', 'MRN', '', 'wcase', '', 'tstp', 0.004, 'type', 1, ...
    'FileNam', file);
for i = 1: 1: header-2
    if regexpi(all{i}, 'PATIENT NAME:')
        hdat.Nam = strtrim(extractAfter(all{i}, 'PATIENT NAME:'));
    elseif regexpi(all{i}, 'PATIENT MRN:')
        hdat.MRN = strtrim(extractAfter(all{i}, 'PATIENT MRN:'));
    elseif regexpi(all{i}, 'WITTCASE NUMBER:')
        hdat.wcase = strtrim(extractAfter(all{i}, 'WITTCASE NUMBER:'));
    elseif regexpi(all{i}, 'TIME RESOLUTION:')
        tmp = strtrim(extractAfter(all{i}, 'TIME RESOLUTION:'));
        tmp = textscan(tmp,'%s',1);
        tmp = tmp{1,1};
        hdat.tstp = str2num(tmp{1})/1000;
    end
end

%% Find columns of RV/LV/ECG traces
lvppos = 0;
ecgpos = 0;

tmp = textscan(all{header+1}, '%s %s %s %s %s %s %s %s',1);
for i = 1 : 2 : 7
    if strcmp(tmp{i},'RV')
        lvppos = i;
    end
    if strcmp(tmp{i},'I')
        if ~strcmp(tmp{i},'III')
            ecgpos = -i;
        else
            ecgpos = i;
        end
    end
end

% If lvppos came back empty, return.
if lvppos == 0 || isempty(lvppos)
    [Pres, dPdt, Rvals, hdat] = deal (0);
    return;
end

%% Process remainder of data as floating point, grabbing ECG indicators
offset = header+1;
Pres = zeros(length(all)-offset,1);
dPdt = Pres;

Pre_trig = [];

for i = header+2 : 1 : length(all)
    % Parse the line into its (possible) components
    DatLin = textscan(all{i}, '%s %s %s %s %s %s %s %s', 1);
    j = i-offset;
    
    % Pressure contains values plus (occasional) ECG markers
    tmp = textscan(char(DatLin{lvppos}), '%f(%c)');
    Pres(j) = tmp{1};
    if ~isempty(tmp{2})
        if strcmp(tmp{2},'S')
            Pre_trig = [Pre_trig j];
            marks = 1;
            if debug, disp(['got an S ' num2str(i)]), end;
        elseif strcmp(tmp{2},'D')
            Pre_trig = [Pre_trig -j];
            marks = 1;
            if debug, disp(['got an D ' num2str(i)]), end;
        end
    end
    
    % dP/dt is just numerical data
    tmp = textscan(char(DatLin{lvppos+1}), '%f');
    if ~isempty(tmp{1});
        dPdt(j) = tmp{1};
    end
    
end

Rvals = abs(Pre_trig(Pre_trig<0));

end
