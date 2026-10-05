clear;
clc;
close all;

set(groot,'defaultFigureVisible','on');

run("config.m");
run("simulate_adc_cube.m");
run("range_fft.m");
run("doppler_fft.m");

