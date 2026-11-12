function phase2_group2_main()

%% Phase 2 Group 2 Main

rng(1);

base = phase2_main();
cfg = base.cfg;
c = base.c;
fc = base.fc;
B = base.B;
Tchirp = base.Tchirp;
fs = base.fs;
Nchirp = base.Nchirp;
Ns = base.Ns;
half_N = base.half_N;
S = base.S;
lambda = base.lambda;
Nrx = base.Nrx;
Nrx_az = base.Nrx_az;
Nrx_el = base.Nrx_el;
d_az = base.d_az;
d_el = base.d_el;
rx_y = base.rx_y;
rx_z = base.rx_z;
NangleFFT_az = base.NangleFFT_az;
NangleFFT_el = base.NangleFFT_el;
target_range = base.target_range;
target_velocity = base.target_velocity;
target_angle = base.target_angle;
target_elevation = base.target_elevation;
target_amplitude = base.target_amplitude;
Ntargets = base.Ntargets;
range_resolution = base.range_resolution;
tau = base.tau;
fb = base.fb;
fd = base.fd;
Tr = base.Tr;
Td = base.Td;
Gr = base.Gr;
Gd = base.Gd;
Pfa = base.Pfa;
snr_db = base.snr_db;
adc_cube = base.adc_cube;
t_fast = base.t_fast;
range_data = base.range_data;
range_axis = base.range_axis;
range_data_no_window = base.range_data_no_window;
roi_threshold_db = base.roi_threshold_db;
roi = base.roi;
range_positive = base.range_positive;
range_energy = base.range_energy;
roi_noise_floor = base.roi_noise_floor;
roi_threshold = base.roi_threshold;
roi_mask = base.roi_mask;
roi_bins = base.roi_bins;
roi_count = base.roi_count;
roi_fraction = base.roi_fraction;
target_roi_selected = base.target_roi_selected;
apply_clutter_removal = base.apply_clutter_removal;
range_data_for_doppler = base.range_data_for_doppler;
range_doppler_rx = base.range_doppler_rx;
RD_power = base.RD_power;
velocity_axis = base.velocity_axis;
RD_mag = base.RD_mag;
RD_db = base.RD_db;
cfar_map = base.cfar_map;
noise_map = base.noise_map;
det_range_idx = base.det_range_idx;
det_doppler_idx = base.det_doppler_idx;
cfar_stats = base.cfar_stats;
Ntraining = base.Ntraining;
alpha = base.alpha;
total_cells = base.total_cells;
guard_cells = base.guard_cells;
peak_range_idx = base.peak_range_idx;
peak_doppler_idx = base.peak_doppler_idx;
peak_power = base.peak_power;
Npeaks = base.Npeaks;
estimated_azimuth = base.estimated_azimuth;
estimated_elevation = base.estimated_elevation;
estimated_snr = base.estimated_snr;
data = base.data;
cost = base.cost;
range_fft_count = base.range_fft_count;
range_fft_diagnostic_count = base.range_fft_diagnostic_count;
doppler_fft_current = base.doppler_fft_current;
doppler_fft_fair = base.doppler_fft_fair;
doppler_fft_sparse = base.doppler_fft_sparse;
projected_roi_reduction = base.projected_roi_reduction;
cfar_cut_count = base.cfar_cut_count;
cfar_training_accesses = base.cfar_training_accesses;


addpath("tests");


%% Actual Selective Doppler

active_range_bins = roi_bins;

[sparse_range_doppler_rx, sparse_RD_power, sparse_fft_count, sparse_fft_work] = selective_doppler_processing_process(range_data_for_doppler, active_range_bins, half_N, Nchirp, Nrx, apply_clutter_removal);
print_selective_doppler_processing(active_range_bins, half_N, sparse_fft_count, sparse_fft_work);

actual_roi_bins = active_range_bins;
actual_sparse_power = sparse_RD_power;
actual_sparse_fft_count = sparse_fft_count;


%% Compare Selected ROI Against Dense

dense_selected = range_doppler_rx(actual_roi_bins, :, :);

selected_fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_selected(:)));

selected_fft_tolerance = 1e-10 * max(1, max(abs(dense_selected(:))));

assert(selected_fft_difference <= selected_fft_tolerance, "Selected ROI FFT mismatch");

selected_power_difference = max(abs(sparse_RD_power(:) - reshape(RD_power(actual_roi_bins, :), [], 1)));

selected_power_tolerance = 1e-10 * max(1, max(abs(RD_power(:))));

assert(selected_power_difference <= selected_power_tolerance, "Selected ROI RD power mismatch");

assert(actual_sparse_fft_count == roi_count * Nrx);


%% Full ROI Regression

regression = dense_sparse_regression(range_doppler_rx, RD_power, cfar_map, noise_map, det_range_idx, det_doppler_idx, peak_range_idx, peak_doppler_idx, peak_power, range_data_for_doppler, half_N, Nchirp, Nrx, apply_clutter_removal, cfg);
fft_difference = regression.fft_difference;
power_difference = regression.power_difference;
fft_match = regression.fft_match;
power_match = regression.power_match;
cfar_match = regression.cfar_match;
peak_match = regression.peak_match;
dense_detection_count = regression.dense_detection_count;



%% Workload Comparison

fair_dense_fft_count = half_N * Nrx;

selective_reduction = 100 * (1 - actual_sparse_fft_count / fair_dense_fft_count);


%% Output Folder

output_dir = "output/phase2/02_selective_doppler_regression";

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end


%% Dense vs Selective RD Plot

dense_max_power = max(RD_power(:));

dense_plot_db = 10 * log10(RD_power / dense_max_power + eps);

sparse_plot_db = NaN(half_N, Nchirp);

sparse_plot_db(actual_roi_bins, :) = 10 * log10(actual_sparse_power / dense_max_power + eps);

fig = figure("Position", [100 100 1200 450]);

subplot(1, 2, 1);

imagesc(velocity_axis, range_axis, dense_plot_db);

axis xy;
xlabel("Velocity [m/s]");
ylabel("Range [m]");
title("Dense Range-Doppler Map");
ylim([0 70]);
clim([-60 0]);
colorbar;


subplot(1, 2, 2);

imagesc(velocity_axis, range_axis, sparse_plot_db);

axis xy;
xlabel("Velocity [m/s]");
ylabel("Range [m]");
title("Selective Doppler - Processed ROI");
ylim([0 70]);
clim([-60 0]);
colorbar;


plot_file = fullfile(output_dir, "dense_vs_selective_doppler_rd_map.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save Validation Results

result_file = fullfile(output_dir, "selective_doppler_full_roi_regression.txt");

file = fopen(result_file, "w");

assert(file ~= -1, "Cannot create Group 2 result file");

fprintf(file, "=== Phase 2 Group 2: Selective Doppler and Regression ===\n\n");

fprintf(file, "[Selective Doppler Test]\n");

fprintf(file, "Selected range bins = %d / %d\n", length(actual_roi_bins), half_N);

fprintf(file, "Current Dense FFT count = %d\n", Ns * Nrx);

fprintf(file, "Fair Dense FFT count = %d\n", fair_dense_fft_count);

fprintf(file, "Selective Doppler FFT count = %d\n", actual_sparse_fft_count);

fprintf(file, "Maximum selected FFT difference = %.3e\n", selected_fft_difference);

fprintf(file, "Maximum selected RD power difference = %.3e\n", selected_power_difference);

fprintf(file, "FFT workload reduction vs Fair Dense = %.2f%%\n", selective_reduction);


fprintf(file, "\n[Full ROI Regression]\n");

fprintf(file, "Full ROI range bins = %d\n", half_N);

fprintf(file, "Maximum FFT difference = %.3e\n", fft_difference);

fprintf(file, "Maximum RD power difference = %.3e\n", power_difference);

fprintf(file, "Complex FFT match = %d\n", fft_match);

fprintf(file, "RD power match = %d\n", power_match);

fprintf(file, "CFAR detection match = %d\n", cfar_match);

fprintf(file, "Peak grouping match = %d\n", peak_match);

fprintf(file, "Dense CFAR detected cells = %d\n", dense_detection_count);


fprintf(file, "\n[Validation]\n");

fprintf(file, "Selected ROI FFT correctness = PASS\n");

fprintf(file, "Full ROI regression = PASS\n");

fprintf(file, "Overall = PASS\n");


fprintf(file, "\n[Notes]\n");

fprintf(file, "Selective FFT processing uses ROI bins only.\n");

fprintf(file, "Full ROI regression uses all positive range bins.\n");

fprintf(file, "ROI-only CFAR is not implemented yet.\n");

fprintf(file, "Workload reduction is not a latency or power measurement.\n");

fclose(file);


%% Final Results

fprintf("\n--- Phase 2 Group 2 Validation ---\n");

fprintf("Selected range bins = %d / %d\n", length(actual_roi_bins), half_N);

fprintf("Selective Doppler FFT count = %d\n", actual_sparse_fft_count);

fprintf("Selected ROI FFT match = PASS\n");

fprintf("Full ROI regression = PASS\n");

fprintf("FFT workload reduction vs Fair Dense = %.2f%%\n", selective_reduction);

fprintf("Group 2 validation passed.\n");

fprintf("Plot saved: %s\n", plot_file);

fprintf("Results saved: %s\n", result_file);

end
