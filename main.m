
clear;
clc;
close all;

addpath(genpath("src"));
set(groot, "defaultFigureVisible", "on");

run("config.m");

run("simulate_adc_cube.m");

run("range_fft.m");

run("static_clutter_removal.m");

run("doppler_fft.m");

run("cfar2d.m");

run("peak_grouping.m");

run("aoa_fft.m");

run("build_pointcloud.m");

run("validate_ground_truth.m");

drawnow;