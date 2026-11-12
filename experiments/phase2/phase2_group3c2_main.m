function phase2_group3c2_main()

%% Phase 2 Group 3C-2 Main

% 함수의 지역 변수 범위가 자동으로 초기화됨
clc;
close all;

assert(isfile("src/config/radar_config.m"), "Run from the project root");

addpath(genpath("src"));
addpath("tests");

rng(1);

cfg = radar_config();
c = cfg.c;
fc = cfg.fc;
B = cfg.B;
Tchirp = cfg.Tchirp;
fs = cfg.fs;
Nchirp = cfg.Nchirp;
target_range = cfg.target_range;
target_velocity = cfg.target_velocity;
target_angle = cfg.target_angle;
target_elevation = cfg.target_elevation;
target_amplitude = cfg.target_amplitude;
Ntargets = cfg.Ntargets;
S = cfg.S;
Ns = cfg.Ns;
range_resolution = cfg.range_resolution;
tau = cfg.tau;
fb = cfg.fb;
lambda = cfg.lambda;
fd = cfg.fd;
Nrx_az = cfg.Nrx_az;
Nrx_el = cfg.Nrx_el;
Nrx = cfg.Nrx;
d_az = cfg.d_az;
d_el = cfg.d_el;
rx_y = cfg.rx_y;
rx_z = cfg.rx_z;
NangleFFT_az = cfg.NangleFFT_az;
NangleFFT_el = cfg.NangleFFT_el;

cfar_pfa_calibration_test();

end
