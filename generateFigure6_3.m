clear all
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

example = 1;
generateFiles = 1;

if generateFiles == 1
    N = 4;
    h = 4;
    curve_length = h;
    save_file = 1;
    Filename1 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvars -except Filename1
    close all
    
    N = 5;
    h = 6;
    curve_length = h;
    save_file = 1;
    Filename2 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvars -except Filename1 Filename2
    close all
    
    N = 6;
    h = 8;
    curve_length = h;
    save_file = 1;
    Filename3 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvar -except Filename1 Filename2 Filename3
    close all
else
    Filename1 = '1_EX4_L';
    Filename2 = '1_EX6_L';
    Filename3 = '1_EX8_L';
end

PlotPictures1;