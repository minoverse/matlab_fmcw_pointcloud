clear;

clc;

close all;


addpath(genpath("src"));


set(groot, "defaultFigureVisible", "on");


run("config.m");


fprintf("Ns = %d samples\n", Ns);

fprintf("Range resolution = %.3f m\n", range_resolution);


for k = 1:Ntargets

    fprintf("\nTarget %d\n", k);

    fprintf("Range = %.2f m\n", target_range(k));

    fprintf("Velocity = %.2f m/s\n", target_velocity(k));

    fprintf("Azimuth = %.2f deg\n", target_angle(k));

    fprintf("Elevation = %.2f deg\n", target_elevation(k));

    fprintf("Beat frequency = %.3f MHz\n", fb(k) / 1e6);

    fprintf("Doppler frequency = %.3f kHz\n", fd(k) / 1e3);

end


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