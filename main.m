clear;
clc;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% Run processing

run("config.m");
run("simulate_adc_cube.m");
run("range_fft.m");
run("doppler_fft.m");

drawnow;


%% Create plots folder

if ~exist("plots", "dir")
    mkdir("plots");
end


%% Save figures after windowing

exportgraphics(figure(1), ...
    "plots/iq_multiple_targets_after_window.png", ...
    "Resolution", 300);

exportgraphics(figure(2), ...
    "plots/range_fft_window_comparison.png", ...
    "Resolution", 300);

exportgraphics(figure(3), ...
    "plots/doppler_target1_after_window.png", ...
    "Resolution", 300);

exportgraphics(figure(4), ...
    "plots/doppler_target2_after_window.png", ...
    "Resolution", 300);

exportgraphics(figure(5), ...
    "plots/rd_map_after_window.png", ...
    "Resolution", 300);