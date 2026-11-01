
%% Phase 2 Group 2 Main

rng(1);

run("phase2_main.m");

addpath("tests");


%% Actual Selective Doppler

active_range_bins = roi_bins;

run("selective_doppler_processing.m");

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

run("dense_sparse_regression.m");


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
