
%% Phase 2 Group 3B Main

clear;
clc;
close all;

assert(isfile("config.m"), "Run from the project root");

addpath(genpath("src"));
addpath("tests");

rng(1);

run("config.m");

run("roi_threshold_miss_rate_test.m");

