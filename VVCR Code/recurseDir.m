% Adam Rauff 
% 6/27/2017
% This function should take a directory name and return all file names in
% the directory and its subdirectories

function [ fileStruct ] = recurseDir(FoldName)

    % make sure input is string and of valid directory
    % ------------------------------------
    
    % get contents of directory name
    dirNam = dir(FoldName);
    
    % initialize structure that holds all files
    fileStruct = dirNam(~[dirNam.isdir]);
    
    % remaining directories
    remDirs = dirNam([dirNam.isdir]);
    for i = 1:length(remDirs)
        if ~strcmp(remDirs(i).name(1),'.')
            dirFiles = recurseDir([FoldName remDirs(i).name filesep]);
            fileStruct = [fileStruct; dirFiles];
        end
    end
end

