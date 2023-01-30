clear all
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

example = 2;
generateFiles = 1;

if generateFiles == 1
    N = 5;
    h = 6;
    curve_length = 15;
    save_file = 1;
    Filename1 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvars -except Filename1
    close all
    
    N = 6;
    h = 8;
    curve_length = 20;
    save_file = 1;
    Filename2 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvars -except Filename1 Filename2
    close all
    
    N = 7;
    h = 10;
    curve_length = 25;
    save_file = 1;
    Filename3 = BFGS_Backtrack(example,N,h,curve_length,save_file);
    clearvar -except Filename1 Filename2 Filename3
    close all
else
    Filename1 = '2_EX15_L';
    Filename2 = '2_EX20_L';
    Filename3 = '2_EX25_L';
end

PlotPictures2;