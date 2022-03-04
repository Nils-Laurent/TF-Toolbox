close all;
addpath('./TF_Toolbox');

addpath('./Validation');

%% downsampling
% MEIGNEN, Sylvain et PHAM, Duong-Hung.
% Retrieval of the modes of multicomponent signals from downsampled short-time fourier transform.
SV_downsamp();

%% ridge detection
SV_exridge();

%% snr dim
SV_snr_dim();

%% sstn
SV_SST();
