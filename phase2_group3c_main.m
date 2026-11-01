
%% Phase 2 Group 3C Main

clear;
clc;
close all;

assert(isfile("config.m"), "Run from the project root");

addpath(genpath("src"));
addpath("tests");

rng(1);

run("config.m");

run("noise_only_monte_carlo_pfa_test.m");
