function [Pres, dPdt, Rvals, file, npath] = ...
    load_calf_p_5_17_17(npath, file, sm)

debug = 0;
marks = 0;

if debug == 1
    disp(['loadp pass ' num2str(sm)]);
end

if isa(file(1),'char')
    filename = strcat(npath, file);   % Get the full filename
    redname = basename(filename, 3);  % Get its basename for reporting
    fd0 = fopen(filename);            % open the file descriptor
    if fd0 < 1
        [Pres, dPdt, Rvals, file, npath] = deal (-1);
        return;
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% AR 5/17
% calf pressure files do not have name or MRN in text file
%     tmp = textscan(fd0, '%s %s %s %s %s %s',1);
%     nam = tmp{3};
%     tmp = textscan(fd0, '%s %s %s',1);
%     mrn = tmp{3};
%     tmp = textscan(fd0, '%s %s %s',1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    tmp = textscan(fd0, '%s %s %s %s %s %s %s',1);
    Time = 0;
    RVPres = 0;
    while ~strcmp (tmp{end},'pressure')
        [tmp, ~] = textscan(fd0, '%s %s %s %s %s %s %s',1);
    end
%     tmp = textscan(fd0, '%s %s %s %s',1);
    for i = 1 : 7
        try
            if strcmp(tmp{i},'Time')
                Time = i;  
            
            elseif strcmp(tmp{i}, 'Millar') && strcmp(tmp{i+1}, 'Pressure') ...
                && strcmp(tmp{i+2}, 'SG')

                MPSG = i;
        
            elseif strcmp(tmp{i}, 'Pressure') && strcmp(tmp{i+1}, 'Systemic')
                RVPres = i-2;
        
            elseif strcmp(tmp{i}, 'Systemic') && strcmp(tmp{i+1}, 'pressure')
                SysPres = i-2;
            end
        catch ME
            Pres = 0; dPdt = 0; Rvals = 0; file = 0; npath = 0;
            return;
        end
    end
    DataALL = textscan(fd0, '%s %s %s %s');  % Read all lines
    
    % check to see "Time" and "RVPres" data are available. some of the text files don't
    % have an "RV" array
    if Time == 0 || RVPres == 0
         
        % set all outputs to 0
        Pres = 0;
        dPdt = 0; 
        Rvals = 0;
        nam = 0;
        mrn = 0;
        file = 0;
        npath = 0;
        marks = 0; 
        
        % exit the function
        disp('Unable to find time array or pressure array');
        return
    end
    
    % extract time array out of cell structure
    T = DataALL{Time};
    Tnumeric = cellfun(@str2double, T);
    
    % extract RV pressure numeric array out of cell structure
    P = DataALL{RVPres};
    Pres = cellfun(@str2double, P);
    
    MP = DataALL{MPSG};
    MPres = cellfun(@str2double, MP);
    
    SP = DataALL{SysPres};
    Spres = cellfun(@str2double, SP);
    
    % computing derivative of pressure by forward difference
    dPdt = diff(Pres)./(Tnumeric(2)-Tnumeric(1));

    % finding the rvals - rough estimates (first round)
    % flip data. used for finding Minima
    DataInv = (-1)*dPdt;
    time_end = Tnumeric(end);
    DerivLen = length(dPdt);
    
    % The peaks must exceed 70 [mmHg/s] and be
    % separated by the number of total seconds of the sample*2+10.
    [pks, pksT] = findpeaks(dPdt, 'MinPeakHeight',70, 'MinPeakDistance',DerivLen/(time_end*2+10));  

    % find peaks of inverted data
    [~, MinIdx] = findpeaks(DataInv, 'MinPeakHeight',70, 'MinPeakDistance',DerivLen/(time_end*2+10));
    
    % check to see if any peaks were found
    if isempty(MinIdx) || isempty(pksT)
        
        Pres = 0;
        dPdt = 0;
        Rvals = 0;
        disp('Signal does not seem to contain apppropriate numbers');
        disp('check load_calf_p file for the pressures and derivative found');
        return
    end
    %finding the local minima values in dP/dt
    Minima = dPdt(MinIdx); 
    
    % data must begin with a dP/dt max, and end with dP/dt min
    if MinIdx(1) < pksT(1)
        MinIdx(1) = [];
        Minima(1) = [];
    end
    if MinIdx(end) < pksT(end)
        pksT(end) = [];
        pks(end) = [];
    end
    
    Rvals = [];
    
    for i = 1:length(MinIdx)
        % find closest pksT that is greater than current min
        tvar = pksT-MinIdx(i);
        tvar(tvar<2) = [];
        tmpMax = min(tvar);
        
        % compute half way between the min and max
        Hway = tmpMax-MinIdx(i);
        
        % convert to time (from index) and add to Rvals vector
        Rvals = [Rvals, Tnumeric(round(MinIdx(i)+Hway))];
    end
    fclose(fd0);
    % disp(['Witt Scientific Pressure File: [...]' redname ' loaded.']);
end
