
%% Phase 2 Group 3A Main

rng(1);

run("phase2_main.m");

addpath("tests");


%% Check Dense Reference

assert(Ntraining == 216);
assert(cfar_cut_count == 25920);


%% Test 1: Actual ROI With CFAR Margin

fprintf("\n=== Group 3A: Selected ROI Test ===\n");

selected_roi_bins = roi_bins;

run("roi_cfar_margin_selection.m");

run("selective_doppler_processing.m");

run("roi_only_cfar.m");

run("roi_cfar_regression.m");


%% Save Selected ROI Results

selected_processing_bins = active_range_bins;

selected_cfar_cuts = roi_cfar_cut_count;

selected_detected_cells = length(roi_det_range_idx);

selected_peak_count = length(roi_peak_range_idx);

selected_det_range_idx = roi_det_range_idx;

selected_det_doppler_idx = roi_det_doppler_idx;

selected_fft_count = sparse_fft_count;

selected_fft_difference = fft_difference;

selected_power_difference = power_difference;

selected_noise_difference = noise_difference;


%% Test 2: Full ROI Regression

fprintf("\n=== Group 3A: Full ROI Test ===\n");

roi_bins = (1:half_N)';

run("roi_cfar_margin_selection.m");

run("selective_doppler_processing.m");

run("roi_only_cfar.m");

run("roi_cfar_regression.m");


%% Full ROI Validation

assert(length(active_range_bins) == half_N);

assert(sparse_fft_count == half_N * Nrx);

assert(roi_cfar_cut_count == cfar_cut_count);

assert(isequal(roi_cfar_local, cfar_map));

full_roi_detected_cells = length(roi_det_range_idx);

assert(full_roi_detected_cells == nnz(cfar_map));


%% Workload Comparison

fair_dense_fft_count = half_N * Nrx;

fft_work_reduction = 100 * (1 - selected_fft_count / fair_dense_fft_count);

cfar_work_reduction = 100 * (1 - selected_cfar_cuts / cfar_cut_count);


%% Output Directory

output_dir = "output/phase2/03_roi_cfar_training_margin";

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end


%% Plot ROI and Margin Coverage

selected_processing_mask = false(half_N, 1);

selected_processing_mask(selected_processing_bins) = true;

fig = figure("Position", [100 100 1100 650]);

subplot(2, 1, 1);

stairs(range_axis, double(selected_processing_mask), "b", "LineWidth", 1.3);

hold on;

stem(range_axis(selected_roi_bins), 1.15 * ones(size(selected_roi_bins)), "r", "filled");

xlabel("Range [m]");
ylabel("Selected");

title("ROI Candidates and CFAR Training Margin");

legend("Processed ROI + Margin", "CFAR CUT Range Bins");

xlim([0 70]);
ylim([0 1.4]);
grid on;


%% Plot ROI CFAR Detections

subplot(2, 1, 2);

imagesc(velocity_axis, range_axis, RD_db);

axis xy;

hold on;

plot(velocity_axis(selected_det_doppler_idx), range_axis(selected_det_range_idx), "go", "MarkerSize", 8, "LineWidth", 1.5);

xlabel("Velocity [m/s]");
ylabel("Range [m]");

title("ROI-Only CFAR Detections (Dense Background for Display)");

ylim([0 70]);
caxis([-60 0]);
colorbar;


%% Save Plot

plot_file = fullfile(output_dir, "roi_cfar_margin_and_detections.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save Validation Summary

result_file = fullfile(output_dir, "roi_cfar_margin_regression.txt");

file = fopen(result_file, "w");

assert(file ~= -1, "Cannot create Group 3A result file");


fprintf(file, "=== Phase 2 Group 3A: ROI-Only CFAR ===\n\n");

fprintf(file, "[Selected ROI]\n");

fprintf(file, "Original ROI bins = %d / %d\n", length(selected_roi_bins), half_N);

fprintf(file, "ROI + CFAR margin bins = %d / %d\n", length(selected_processing_bins), half_N);

fprintf(file, "CFAR range margin = %d bins per side\n", Tr + Gr);

fprintf(file, "Fair Dense Doppler FFT count = %d\n", fair_dense_fft_count);

fprintf(file, "Sparse FFT count with margin = %d\n", selected_fft_count);

fprintf(file, "Doppler FFT workload reduction vs Fair Dense = %.2f%%\n", fft_work_reduction);


fprintf(file, "\n[ROI-Only CFAR]\n");

fprintf(file, "Training cells = %d\n", Ntraining);

fprintf(file, "Dense CFAR CUTs = %d\n", cfar_cut_count);

fprintf(file, "ROI CFAR CUTs = %d\n", selected_cfar_cuts);

fprintf(file, "ROI CFAR workload reduction = %.2f%%\n", cfar_work_reduction);

fprintf(file, "ROI detected cells = %d\n", selected_detected_cells);

fprintf(file, "ROI grouped peaks = %d\n", selected_peak_count);


fprintf(file, "\n[Selected ROI Regression]\n");

fprintf(file, "Maximum FFT difference = %.3e\n", selected_fft_difference);

fprintf(file, "Maximum RD power difference = %.3e\n", selected_power_difference);

fprintf(file, "Maximum CFAR noise difference = %.3e\n", selected_noise_difference);

fprintf(file, "FFT / RD / CFAR / Peak regression = PASS\n");


fprintf(file, "\n[Full ROI Regression]\n");

fprintf(file, "Full ROI bins = %d\n", half_N);

fprintf(file, "Full ROI CFAR CUTs = %d\n", roi_cfar_cut_count);

fprintf(file, "Full ROI detected cells = %d\n", full_roi_detected_cells);

fprintf(file, "Full ROI regression = PASS\n");


fprintf(file, "\n[Notes]\n");

fprintf(file, "Margin FFT workload is included in the sparse cost.\n");

fprintf(file, "CFAR detections are evaluated at ROI CUTs only.\n");

fprintf(file, "Full ROI regression verifies equality with Dense.\n");

fprintf(file, "Dense RD background is used only for plotting.\n");

fprintf(file, "No real latency or power reduction is claimed.\n");


fprintf(file, "\nOverall = PASS\n");

fclose(file);


%% Restore Original ROI Selection

roi_bins = selected_roi_bins;


%% Final Results

fprintf("\n--- Phase 2 Group 3A Validation ---\n");

fprintf("Original ROI bins = %d\n", length(selected_roi_bins));

fprintf("ROI + margin bins = %d\n", length(selected_processing_bins));

fprintf("Sparse Doppler FFTs including margin = %d\n", selected_fft_count);

fprintf("ROI CFAR CUTs = %d / %d\n", selected_cfar_cuts, cfar_cut_count);

fprintf("FFT workload reduction vs Fair Dense = %.2f%%\n", fft_work_reduction);

fprintf("CFAR CUT workload reduction = %.2f%%\n", cfar_work_reduction);

fprintf("Selected ROI regression = PASS\n");

fprintf("Full ROI regression = PASS\n");

fprintf("Group 3A validation passed.\n");

fprintf("Plot saved: %s\n", plot_file);

fprintf("Results saved: %s\n", result_file);

