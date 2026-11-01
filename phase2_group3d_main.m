
%% Phase 2 Group 3D Main

clear;
clc;
close all;

assert(isfile("config.m"), "Run from the project root");

addpath(genpath("src"));
addpath("tests");

assert(isfile("src/optimization/roi_threshold_gate.m"));
assert(isfile("src/optimization/selective_doppler_processing.m"));

rng(2);

run("config.m");

run("cfar_calibrated_detection_test.m");
