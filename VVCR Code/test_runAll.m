function testAll
% Based on runAll, Adam Rauff 10.1.2016
% This code performs multiple calls to calc_vvcr (old VVCR_MULTIH) while varying
% one of the landmarks to determine how Pmax varies with this landmark.

% Prepare the workspace
close all;
clc;

for isopos = 1: 1: 4
% This is the parameter to vary: 1-Ps, 2-Pe, 3-Ns, 4-Ne.
%isopos = 1;

% bring up UI that allows user to select a folder (directory)
%Fold_name = uigetdir('','Select folder containing txt files with pressure data');
%if ~Fold_name
%    disp('runAll: uigetdir closed.');
%    return;
%end
Fold_name = 'C:\Users\hunterk\Desktop\Students\Akshay\20181002_TedfordData\FakeWITT';

% Output Command Window info to diary file
diary('RV_stdout.txt');
mydate = datestr(now,'ddd mmm dd HH:MM:SS YYYY');
disp(['VVCR Analysis ' mydate]);
disp(['              ' Fold_name]);

% attach a forward slash to the end of Fold_name so folder address ends
% with /
Fold_name = [Fold_name filesep];

% list all filenames in the directory of the chosen folder and its
% subdirectories
top = recurseDir(Fold_name);

% keep track of how many pressure files were analyzed
fileCount = 0;

% keep track of skipped files
FileSkpTrck = [];

% keep track of analyzed files
NewAnalyzedTxt = [];

% create variable to hold name of csv file
% In future: should add an input box when script runs so user can interactively
% assign csv file name every time script runs
csvName = ['RV_isopos' num2str(isopos) '.csv'];

agrName = ['RV_isopos' num2str(isopos) '.agr'];
plt_fh = fopen(agrName, 'w');
if ~plt_fh
    error('agr file failed to open, unknown error');
end

% Here would be a good opportunity to open the csv file if it exists, and check
% if top_name already exists in the csv file. If it already exists in the csv
% file, the analysis was already run on it, and it should be skipped check to if
% csvName file exists
filechk = which(csvName);

if ~isempty(filechk)
    fd0 = fopen(csvName);

    % read all columns of csv file (11 columns)
    csvDat = textscan(fd0, '%s %s %s %s %s %s %s %s %s %s %s', ...
        'HeaderLines',1, 'Delimiter',',');
    
    OldAnalyzedTxt = csvDat{1,3}; % this will be a Ax1 cell with A = number of text files analyzed
    
    fclose(fd0);
else
    OldAnalyzedTxt = [];
end

for i = 1:length(top)
    % get individual file (ith file)
    top_name = top(i).name;
    Fold_name = [top(i).folder filesep];
    
    % -----------------------------------------------------------------    
    % this piece of code is only added for the analysis of the text files that
    % were recieved from Dr. Hunter on 9/30/16. None of these files had a .txt
    % file extension, and therefore never satisfied the regexp if statement
    % below.
    [pathstr, name, ext] = fileparts(strcat(Fold_name,'/',top_name));
    
    if ~isempty(strfind('H', top_name)) && ~strcmpi(ext, '.txt')
       top_name = strcat(top_name, '.txt'); 
    end
    % -------------------------------------------------------------------
    
    % check to see if .txt file, then run VVCR analysis
    if regexp (top_name, 'txt$')
        
        TXT_FLAG = false; % this flag turns true after VVCR_* succesfully ran
        % there could be a scenario of a .txt file not being a pressure
        % data file, a corrupted file, empty, etc.
        
        SKIP_FLAG = false; % this flag turn true when the current (ith) file has already been analyzed
        
        if ~isempty(OldAnalyzedTxt)
            for j = 1: length(OldAnalyzedTxt)
                % if the file has been found in the list of analyzed files,
                % skip this iteration of i and continue to the next i
                if strcmp(OldAnalyzedTxt{j},top_name)
                    disp('File has already been analyzed in previous session!');
                    disp(['name of file: ', top_name]);
                    SKIP_FLAG = true;
                end
            end
        end
        
        if ~isempty(NewAnalyzedTxt)
            for j = 1:size(NewAnalyzedTxt,1)
                if strcmp(NewAnalyzedTxt(j,:), top_name)
                    disp('File has already been analyzed THIS session!');
                    disp(['name of file: ', top_name]);
                    SKIP_FLAG = true;
                end
            end
        end
        
        if SKIP_FLAG == true
            % go to next iteration of i
            continue
        end

        % Ignore silently output & info files stored in data directory.
        if ~isempty(regexp(top_name, 'RV_stdout')) || ...
                ~isempty(regexp(top_name, 'README'))
            continue
        end

        % print current file being analyzed
        disp(['Current File: ', top_name]);
        
        disp('TIMING2 ---------------------------');
        % add try catch here. In case txt file is not pressure data, or what
        % not. In the catch, record the name of the text file somewhere, so if
        % it is a pressure file that could not be opened, it is remarked Did
        % this within the load_calf file, it is a good idea!! Thanks Adam!
        [Res0, Pat0] = test_calc_vvcr(Fold_name, top_name, 0, isopos);
        [Res1, Pat1] = test_calc_vvcr(Fold_name, top_name, 1, isopos);
        [Res2, Pat2] = test_calc_vvcr(Fold_name, top_name, 2, isopos);
        [Res3, Pat3] = test_calc_vvcr(Fold_name, top_name, 3, isopos);
        [Res4, Pat4] = test_calc_vvcr(Fold_name, top_name, 4, isopos);

        if ~isstruct(Res2)
            % check to see if outputs are false or true
            % EXIT
            if Res1 == false
                % exit for loop
                % keep in mind when the exit button is pressed, the current
                % patient, i, will not be evaluated
                if fileCount >= 1
                    disp(['The last file analyzed was: ', top(i-1).name]);
                else
                    disp('No File was Analyzed in this session');
                end

                break

            % Discard Patient
            elseif Res2 == true
                % proceed to next iteration of i (next patient)

                disp(['File ',top(i).name, ' is skipped']);

                % write to csv

                filechk = which(csvName);

                % if file does not exist, create one & write header line.
                if isempty(filechk)
                    fd0 = fopen(csvName, 'w');
                    fprintf(fd0, header);
                else
                    % if Pat.FileNam exists, append to it
                    fd0 = fopen(csvName, 'a');

                    if i == 1 || fileCount == 0 
                        fprintf(fd0, ' ,\n');
                    end
                end

                myStr = ['This, Pat.FileNam, ',top(i).name, ...
                    ', was, skipped, becuase, user, determined, ' ...
                    'something, was, wrong'];
                fprintf(fd0,[myStr, '\n']);

                % close file
                fclose(fd0);

                % remark iteration of i that was skipped
                FileSkpTrck = i;

                % skip to next iteration of i
                continue 
            end
        end
        
        % At this point we know that the return is valid (is a structure).
        % We build the file header based on the returned fields in Res,
        % ignoring (for now) the fields that are structures (these are the
        % total results of each analysis.
        isoposFields = {'iso1Ps_Mean', 'iso2Pe_Mean', 'iso3Ns_Mean', ...
            'iso4Ne_Mean'};
        isoposDFields = {'iso1dPs_Mean', 'iso2dPe_Mean', 'iso3dNs_Mean', ...
            'iso4dNe_Mean'};
        deltP = 'DELTA_Pmax_Mean';
        deltV = ['DELTA' isoposFields{isopos}];
        
        resh = {'C_KindNorm_Pmax_Mean', isoposFields{isopos}, ...
            isoposDFields{isopos}, deltP, deltV, 'Deriv'};

        Res1.(deltP) = 0;
        Res1.(deltV) = 0;
        Res1.Deriv = 0;
        
        header = ['file,'];
        for i = 1: 1: length(resh)
            field = resh{i};
            if ~isstruct(Res1.(field))
                % At this point we know it's going to be output, and things have
                % been alphabetized as desired, but we perform transformations
                % on the field names to make them a little more human readable.
                field = regexprep(field, '^[A-Z]_', '');
                field = regexprep(field, '[1-9]', '');
                field = regexprep(field, '_Vcorr', ' Vcorr');
                field = regexprep(field, '_nFit', ' nFit');
                if length(regexp(field, '_')) > 1
                    field = regexprep(field, '_', ' ','once');
                end
                header = [header field ','];
            end
        end
        header = [header '\r\n'];
        
        % check Pat.Nam and Pat.MRN, to make sure they do not give valid
        % names and mrn. The analyzed data must be anonymized!!!

        TXT_FLAG = true; % text file was readable by loadp
        fileCount = fileCount + 1;

        % record the indice of the file associated with each fileCount
        FileCountIter(fileCount) = i;

        % record the name of the file. this is used to make sure 2 of
        % the same files aren't analyzed
        if isempty(NewAnalyzedTxt)
            NewAnalyzedTxt = top_name;
        else
            NewAnalyzedTxt = [NewAnalyzedTxt; top_name];
        end

        if TXT_FLAG == true
            % check to if csvName file exists
            filechk = which(csvName);

            % if file does not exist, create one & write header line.
            if isempty(filechk)
                fd0 = fopen(csvName, 'w');
                fprintf(fd0, header);
            else
                % if file exists, append to it
                fd0 = fopen(csvName, 'a');
                
                % After a subsequent session begins, the program appended data
                % to the last filled cell. This statement makes sure to start a
                % new line before printing the data, to ensure the appended data
                % is recorded separatly, and does not corrupt the data of the
                % last cell
                if i == 1 || fileCount == 1
                    fprintf(fd0, ' ,\r\n');
                end
            end

            % Name, MRN, Filename
            NamNoComma = regexprep (Pat2.Nam, ',', '');
            fprintf(fd0, '%s,', Pat1.FileNam);

            % print remainder based entirely on field names avoiding (for
            % now) returned structures.
            Res = Res0;
            fprintf(plt_fh, '@type xy\r\n');
            for j = 0: 1: 4
                Pre = Res;
                eval(['Res = Res' num2str(j), ';']);

                Res.(deltP) = Res.C_KindNorm_Pmax_Mean - Pre.C_KindNorm_Pmax_Mean;
                Res.(deltV) = Res.(isoposFields{isopos}) - Pre.(isoposFields{isopos});
                if j > 0
                    fprintf(fd0, ' ,');
                    Res.Deriv = Res.(deltP)/Res.(deltV);
                else
                    Res.Deriv = 0;
                end
                
                for i = 1: 1: length(resh)
                    field = resh{i};
                    if ~isstruct(Res.(field))
                        value = Res.(field);
                        if round(value) == value
                            fprintf(fd0, '%4i, ', value);
                        else
                            fprintf(fd0, '%10.6f, ', value);
                        end
                    end
                end

                fprintf(fd0, '\r\n');

                fprintf(plt_fh, '%12.8f %12.8f\r\n', Res.(deltV), Res.(deltP));
            end

            fprintf(plt_fh, '&\r\n');

            % close Pat.FileNam
            fclose(fd0);
            
            % print to command window how many files were analyzed
            disp(['FileCount: ', num2str(fileCount)]);
            
        end   
    end
end

end

fclose all;
diary off;

end
