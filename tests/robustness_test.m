clear;
clc;
close all;

project_root = fileparts(fileparts(mfilename("fullpath")));
cd(project_root);

addpath(genpath("src"));
addpath("tests");

%% Choose test

% 1 = Weak Target
% 2 = Close Targets
% 3 = Static Clutter
% 4 = RX Gain / Phase Mismatch
% 5 = Near-Zero Velocity
% 6 = Elevation Check

test_case = 6;


%% Near-zero setting

near_zero_notch = true;


run("config.m");


apply_clutter_removal = false;


%% Test 1 - Weak Target

if test_case == 1

    target_amplitude = [1.0, 0.2];

end


%% Test 2 - Close Targets

if test_case == 2

    target_range = [30, 30.45];

    target_velocity = [10, 10];

    target_angle = [15, -20];

    target_elevation = [10, -8];

    target_amplitude = [1.0, 0.7];

end


%% Test 3 - Static Clutter

if test_case == 3

    target_range = [30, 55, 20];

    target_velocity = [10, -5, 0];

    target_angle = [15, -20, 0];

    target_elevation = [10, -8, 0];

    target_amplitude = [1.0, 0.7, 2.0];

    apply_clutter_removal = true;

end


%% Test 5 - Near-Zero Velocity

if test_case == 5

    target_range = [30, 55, 40];

    target_velocity = [10, -5, 0.6];

    target_angle = [15, -20, 0];

    target_elevation = [10, -8, 5];

    target_amplitude = [1.0, 0.7, 1.0];

    apply_clutter_removal = near_zero_notch;

end


%% Test 6 - Elevation Check

if test_case == 6

    target_elevation = [20, -15];

end


%% Recalculate target parameters

Ntargets = length(target_range);

tau = 2 * target_range / c;

fb = S .* tau;

fd = 2 * target_velocity / lambda;


%% Generate ADC

run("simulate_adc_cube.m");


%% Test 4 - RX Gain / Phase Mismatch

if test_case == 4

    rx_gain = [1.0, 0.9, 1.1, 0.85];

    rx_phase = [0, 8, -5, 12];


    for rx = 1:Nrx

        phase_error = exp(1j * rx_phase(rx) * pi / 180);

        adc_cube(:, rx, :) = adc_cube(:, rx, :) * rx_gain(rx) * phase_error;

    end

end


%% Processing pipeline

run("range_fft.m");

run("static_clutter_removal.m");

run("doppler_fft.m");

run("cfar2d.m");

run("peak_grouping.m");

run("aoa_fft.m");

run("build_pointcloud.m");

run("validate_ground_truth.m");


%% Save

if ~exist("plots/robustness", "dir")

    mkdir("plots/robustness");

end


if ~exist("output", "dir")

    mkdir("output");

end


if test_case == 5

    if near_zero_notch

        filename = "plots/robustness/test_5_notch_on_pointcloud.png";

        csv_name = "output/robustness_test_5_notch_on.csv";

    else

        filename = "plots/robustness/test_5_notch_off_pointcloud.png";

        csv_name = "output/robustness_test_5_notch_off.csv";

    end

else

    filename = "plots/robustness/test_" + test_case + "_pointcloud.png";

    csv_name = "output/robustness_test_" + test_case + ".csv";

end


exportgraphics(gcf, filename, "Resolution", 300);

writetable(pointcloud, csv_name);