
%% Phase 2 Group 3C-2 Main

clear;
clc;
close all;

assert(isfile("config.m"), "Run from the project root");

addpath(genpath("src"));
addpath("tests");

rng(1);

run("config.m");

run("cfar_pfa_calibration_test.m");
