% Run VVCR analysis on a list of text files
% Adam Rauff 10.1.2016

% Prepare the workspace
clear all;
close all;
clc;

% bring up UI that allows user to select a folder (directory)
Fold_name = uigetdir('','Select folder containing txt files with pressure data');

% attach a forward slash to the end of Fold_name so folder address ends
% with /
Fold_name = [Fold_name, '/'];

% list all the stuff in the directory of the chosen folder
top  = dir(Fold_name);

% keep track of how many pressure files were analyzed
fileCount = 0;

% keep track of skipped files
FileSkpTrck = [];

for i = 1:length(top)
    
    % get individual file (ith file)
    top_name = top(i).name;
    
    % check to see if .txt file, then run VVCR analysis
    if regexp (top_name, 'txt$')
        
        TXT_FLAG = false; % this flag turns true after VVCR_FINAL succesfully ran
        % there could be a scenario of a .txt file not being a pressure
        % data file, a corrupted file, empty, etc.
        
        % add try catch here. In case txt file is not pressure data, or
        % what not. In the catch, record the name of the text file
        % somewhere, so if it is a pressure file that could not be opened,
        % it is remarked
%         try

            % Here would be a good opportunity to open the csv file if it
            % exists, and check if top_name already exists in the csv file.
            % If it already exists in the csv file, the analysis was
            % already run on it, and it should be skipped
            % check to if RV_data.csv file exists
            filechk = which('RV_data.csv');
            
%             if ~isempty(filechk)
%                 fd0 = fopen('RV_data.csv');
%                 fclose(fd0);
%             end
            [AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH, Pnam, Pmrn, file] = VVCR_FINAL_10_3(Fold_name,top_name);
            
            % check to see if outputs are false or true
            % EXIT
            if AVG_Pes == false
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
            elseif AVG_Pes == true
                % proceed to next iteration of i (next patient)
                
                disp(['File ',top(i).name, ' is skipped']);
                
                % write to csv
                
                filechk = which('RV_data.csv');

                % if file does not exist, creat one, and write all the headers to it
                if isempty(filechk)
                    fd0 = fopen('RV_data.csv', 'w');
                    fprintf(fd0, 'Pnam, Pmrn, file, AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH\n');
                else
                    % if file exists, append to it
                    fd0 = fopen('RV_data.csv', 'a');
                end
                
                myStr = ['This, file, ',top(i).name, ', was, skipped'];
                fprintf(fd0,[myStr, '\n']);
                
                % close file
                fclose(fd0);
                
                % remark iteration of i that was skipped
                FileSkpTrck = i;
                
                % skip to next iteration of i
                continue 
            end
            % ----------------------------------------------------------
            % check Pnam and Pmrn, to make sure they do not give valid
            % names and mrn. The analyzed data must be anonymized!!!
            % ----------------------------------------------------------
            
            TXT_FLAG = true; % text file was readable by loadp
            fileCount = fileCount + 1;
            
            % record the indice of the file associated with each fileCount
            FileCountIter(fileCount) = i;
%         catch
%             disp('file could not be analyzed');
%             disp(['name of file: ',top_name]);
            
%         end
        
        if TXT_FLAG == true
            % check to if RV_data.csv file exists
            filechk = which('RV_data.csv');

            % if file does not exist, creat one, and write all the headers to it
            if isempty(filechk)
                fd0 = fopen('RV_data.csv', 'w');
                fprintf(fd0, 'Pnam, Pmrn, file, AVG_Pes, AVG_Pmax, VVCR_UT, VVCR_KH\n');
            else
                % if file exists, append to it
                fd0 = fopen('RV_data.csv', 'a');
                fprintf(fd0, ' ,\n');
            end
            
            % After a file is skipped, the program appended data to the
            % last filled cell. This statement makes sure to start a new
            % line before printing the data, to ensure the appended data is
            % recorded separatly, and does not corrupt the data of the last
            % cell
%             if ~isempty(FileSkpTrck)
%                 if FileSkpTrck == i-1 
%                     fprintf(fd0, ' ,\n');  
%                 end
%             end
            
            % Name, MRN, Filename
            fprintf(fd0, '%s, %s, %s,', Pnam{1}, Pmrn{1}, file);
                    
            % average Pes and Pmax - end systolic pressure (Pes), and
            % maximum pressure (Pmax) obtained from equation of Naeiji
            % (Mean + 2*amp)
            fprintf(fd0, '%10.6f, %10.6f,' , AVG_Pes, AVG_Pmax);
            
            % Finally, VVCR - ventricular vascular coupling ratio
            % UT - Dr. Uyen Troung
            % KH - Dr. Kendall Hunter
            % VVCR_KH = 1/VVCR_UT, they are reciprocals
            fprintf(fd0, '%8.6f, %8.6f\n', VVCR_UT, VVCR_KH);
            
            % close file
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