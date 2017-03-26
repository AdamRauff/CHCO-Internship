% 3/14/2016
% XML parser for TomTec Output files 
% Adam Rauff & KSH :-)

% add anaconda path to matlab search path

% path('/Users/adamrauff/anaconda/bin',path)   % OSX
addpath('C:\Program Files\Anaconda3\Scripts')  % Windows 10 x64


% run python script
% system allows matlab to access bash shell like enironment
% provide - 'source anaconda-path/activate env-name;python scriptname args'

%system('source /Users/adamrauff/anaconda/bin/activate OpenCV;python RunPXML.py TomTec\ Data/')  % OSX
system('activate.bat OpenCV & python RunPXML.py "TomTec Data/"')                                 % Windows 10 x64