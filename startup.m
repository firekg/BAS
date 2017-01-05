% a general startup script:
disp(['Adding all subfolders to path...']);
me = mfilename;
mydir = which(me);
mydir = mydir(1:end-2-numel(me));
addpath(genpath(mydir(1:end-1)))
clear me mydir
disp(['Adding all subfolders to path... Done.']);
% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup'), PsychStartup; end;

