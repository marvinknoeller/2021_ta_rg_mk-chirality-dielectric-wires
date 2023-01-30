clear all
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

example = 1;
generateFiles = 1;

if generateFiles == 1
    N = 5;
    h = 6;
    curve_length = h;
    save_file = 1;
    Filename = BFGS_Backtrack(example,N,h,curve_length,save_file);
else
    Filename = '1_EX6_L';
end


%% in case that the .mat file has already been generated
clearvars -except Filename
load(strcat(Filename,'.mat'));
Plot_steps_vs_measure_1;
% epsilonscan first generates, then plots
clearvars -except Filename
load(strcat(Filename,'.mat'));
EpsilonScan1;
% wavenumberscan first generates, then plots
clearvars -except Filename
load(strcat(Filename,'.mat'));
WavenumberScan1;