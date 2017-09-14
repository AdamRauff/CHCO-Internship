% Run VVCR analysis on a list of text files
% Adam Rauff 10.1.2016

% Prepare the workspace
clear all;
close all;
clc;

% bring up UI that allows user to select a folder (directory)
Fold_name = uigetdir('','Select folder containing txt files with pressure data');
if ~Fold_name
    disp('runAll: uigetdir closed.');
    return;
end

% attach a forward slash to the end of Fold_name so folder address ends
% with /
Fold_name = [Fold_name, '/'];

% list all filenames in the directory of the chosen folder and its
% subdirectories
top  = recurseDir(Fold_name);

% keep track of how many pressure files were analyzed
fileCount = 0;

% keep track of skipped files
FileSkpTrck = [];

% keep track of analyzed files
NewAnalyzedTxt = [];

% create variable to hold name of csv file
% In future: should add an input box when script runs so user can
% interactively assign csv file name every time script runs
csvName = 'RV_data_03.csv';

% Here would be a good opportunity to open the csv file if it
% exists, and check if top_name already exists in the csv file.
% If it already exists in the csv file, the analysis was
% already run on it, and it should be skipped
% check to if RV_data.csv file exists
filechk = which(csvName);

if ~isempty(filechk)
    fd0 = fopen(csvName);

    % read all columns of csv file (11 columns)
    csvDat = textscan(fd0, '%s %s %s %s %s %s %s %s %s %s %s', 'HeaderLines',1, 'Delimiter',',');
    
    OldAnalyzedTxt = csvDat{1,3}; % this will be a Ax1 cell with A = number of text files analyzed
    
    fclose(fd0);
else
    OldAnalyzedTxt = [];
end

for i = 1:length(top)
    % get individual file (ith file)
    top_name = top(i).name;
    
    % -----------------------------------------------------------------
    % this piece of code is only added for the analysis of the text files
    % that were recieved from Dr. Hunter on 9/30/16. None of these files
    % had a .txt file extension, and therefore never satisfied the regexp
    % if statement below.
    [pathstr, name, ext] = fileparts(strcat(Fold_name,'/',top_name));
    
    if ~isempty(strfind('H', top_name)) && ~strcmpi(ext, '.txt')
       top_name = strcat(top_name, '.txt'); 
    end
    % -------------------------------------------------------------------
    
    % check to see if .txt file, then run VVCR analysis
    if regexp (top_name, 'txt$')
        
        TXT_FLAG = false; % this flag turns true after VVCR_*.m succesfully ran
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
        
        % print current file being analyzed
        disp(['Current File: ', top_name]);
        
        % add try catch here. In case txt file is not pressure data, or
        % what not. In the catch, record the name of the text file
        % somewhere, so if it is a pressure file that could not be opened,
        % it is remarked
%         try
            [Res, Pat] = VVCR_MULTIH_08_09_17(Fold_name,top_name);
            
            if ~isstruct(Res)
                % check to see if outputs are false or true
                % EXIT
                if Res == false
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
                elseif Res == true
                    % proceed to next iteration of i (next patient)

                    disp(['File ',top(i).name, ' is skipped']);

                    % write to csv

                    filechk = which(csvName);

                    % if file does not exist, creat one, and write all the 
                    % headers to it
                    if isempty(filechk)
                        fd0 = fopen(csvName, 'w');
                        fprintf(fd0, ['Pnam, Pmrn, file, AVG_Pes, ' ...
                            'AVG_Pmax, VVCR_UT, VVCR_KH, Num_Peaks, ' ...
                            'STD_Pes, STD_PMX, TotNumWaves\n']);
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
            % ----------------------------------------------------------
            % check Pat.Nam and Pat.MRN, to make sure they do not give
            % valid names and mrn. The analyzed data must be anonymized!!!
            % ----------------------------------------------------------
            
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
            
%         catch
%             disp('file could not be analyzed');
%             disp(['name of Pat.FileNam: ',top_name]);
            
%         end
        
        if TXT_FLAG == true
            % check to if RV_data.csv Pat.FileNam exists
            filechk = which(csvName);

            % if file does not exist, creat one, and write all the headers to it
            if isempty(filechk)
                fd0 = fopen(csvName, 'w');
                fprintf(fd0, 'Pnam, Pmrn, file, AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Num_Peaks, STD_Pes, STD_PMX, TotNumWaves\n');
            else
                % if file exists, append to it
                fd0 = fopen(csvName, 'a');
                
                % After a subsequent session begins, the program appended data to the
                % last filled cell. This statement makes sure to start a new
                % line before printing the data, to ensure the appended data is
                % recorded separatly, and does not corrupt the data of the last
                % cell
                if i == 1 || fileCount == 1
                    fprintf(fd0, ' ,\n');
                end
            end
            
            % Name, MRN, Filename
            fprintf(fd0, '%s, %s, %s,', Pat.Nam{1}, Pat.MRN{1}, Pat.FileNam);
                    
            % average Pes and Pmax - end systolic pressure (Pes), and
            % maximum pressure (Pmax) obtained from equation of Naeiji
            % (Mean + 2*amp)
            fprintf(fd0, '%10.6f, %10.6f,' , Res.Pes_Mean, Res.Pmax_Mean);
            
            % print VVCR - ventricular vascular coupling ratio
            % UT - Dr. Uyen Troung
            % KH - Dr. Kendall Hunter
            % Res.VVCRnorm_Mean = 1/Res.VVCRinv_Mean, they are reciprocals
            fprintf(fd0, '%8.6f, %8.6f,', Res.VVCRinv_Mean, Res.VVCRnorm_Mean);
            
            % print the standard deviations of Pmax and Pes, the number
            % of analyzed peaks, and the total number of waves
            fprintf(fd0, '%i, %10.6f, %10.6f, %i\n', Res.numPeaks, Res.Pes_StD, Res.Pmax_StD, Res.TotNumWaves);
            
            % close Pat.FileNam
            fclose(fd0);
            
            % print to command window how many files were analyzed
            disp(['FileCount: ', num2str(fileCount)]);
%         else
%             % check to if RV_data.csv file exists
%             filechk = which('RV_data.csv');
%             
%             % if file does not exist, creat one, and write all the headers to it
%             if isempty(filechk)
%                 fd0 = fopen('RV_data.csv', 'w');
%                 fprintf(fd0,'%s did not load\n', top_name);
%                 fclose(fd0);
%             else
%                 % if file exists, append to it
%                 fd0 = fopen('RV_data.csv', 'a');
%                 fprintf(fd0,'%s did not load\n', top_name);
%                 fclose(fd0);
%             end
            
        end   
    end
end
