%% clear everything
close all
clear
clc

%% read input data
% images = {'images/lena.bmp', 'images/peppers.bmp', 'images/boat.bmp'};
images = {'images/lena.bmp'};
cI = readImages(images);
cI = cI./255;

%% parameters
% number of non-zeros per column in the sparse factor S of the dictionary
s = 3;
% the number of G transforms in the orthonormal factor U of the dictionary
g = 170;
% the sparsity
k0 = 4;
% size of the dictionary D (number of atoms)
m = 128;

%% calls
% initialization with DCT
[U_dct, S_dct, supportS, X_dct, tus_dct, err_dct] = dct_f_dla(cI, k0, m, s);
% F-DLA
[U_fdla, S_fdla, X_fdla, positions_fdla, values_fdla, tus_fdla, err_fdla] = f_dla(cI, k0, g, m, U_dct, S_dct, supportS, X_dct);
