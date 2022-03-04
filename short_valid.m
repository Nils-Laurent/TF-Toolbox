close all;
addpath('./TF_Toolbox');

addpath('./Validation');

%% downsampling
SV_downsamp();

%% ridge detection
SV_exridge();

%% snr dim
SV_snr_dim();

%% sstn
SV_SST();

%% add noise
SV_addNoise();
