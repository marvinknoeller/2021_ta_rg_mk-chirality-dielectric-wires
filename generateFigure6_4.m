clear all
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

example = 2;
generateFiles = 1;

if generateFiles == 1
    N = 6;
    h = 8;
    curve_length = 20;
    save_file = 1;
    Filename = BFGS_Backtrack(example,N,h,curve_length,save_file);
else
    Filename = '2_EX20_L';
end


load(strcat(Filename,'.mat'));
Plot_Iterates_Example2;