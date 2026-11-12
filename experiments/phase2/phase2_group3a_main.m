function phase2_group3a_main()

%% Phase 2 Group 3A Main

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


%% Check Dense Reference

assert(Ntraining == 216);
assert(cfar_cut_count == 25920);


%% Test 1: Actual ROI With CFAR Margin

fprintf("\n=== Group 3A: Selected ROI Test ===\n");

selected_roi_bins = roi_bins;

roi_margin = roi_cfar_margin_selection_process(roi_bins, half_N, Tr, Gr);
roi_range_margin = roi_margin.roi_range_margin;
roi_margin_mask = roi_margin.roi_margin_mask;
active_range_bins = roi_margin.active_range_bins;
roi_processing_count = roi_margin.roi_processing_count;
roi_extra_margin_bins = roi_margin.roi_extra_margin_bins;
print_roi_cfar_margin_selection(roi_margin, half_N, Nrx);

[sparse_range_doppler_rx, sparse_RD_power, sparse_fft_count, sparse_fft_work] = selective_doppler_processing_process(range_data_for_doppler, active_range_bins, half_N, Nchirp, Nrx, apply_clutter_removal);
print_selective_doppler_processing(active_range_bins, half_N, sparse_fft_count, sparse_fft_work);

radar_roi = struct("half_N",half_N,"Nchirp",Nchirp,"Tr",Tr,"Td",Td,"Gr",Gr,"Gd",Gd,"Pfa",Pfa,"Ntraining",Ntraining);
[roi_cfar_local, roi_noise_local, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, roi_stats, roi_det_range_idx, roi_det_doppler_idx, range_lookup] = roi_only_cfar_process(roi_bins, active_range_bins, sparse_RD_power, radar_roi);
print_roi_only_cfar(roi_stats);

regression = roi_cfar_regression(range_doppler_rx, active_range_bins, sparse_range_doppler_rx, RD_power, sparse_RD_power, cfar_map, roi_bins, roi_cfar_local, noise_map, roi_noise_local, peak_range_idx, peak_doppler_idx, roi_peak_range_idx, roi_peak_doppler_idx);
fft_difference = regression.fft_difference;
power_difference = regression.power_difference;
noise_difference = regression.noise_difference;
fft_match = regression.fft_match;
power_match = regression.power_match;
cfar_match = regression.cfar_match;
noise_match = regression.noise_match;
peak_match = regression.peak_match;



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

roi_margin = roi_cfar_margin_selection_process(roi_bins, half_N, Tr, Gr);
roi_range_margin = roi_margin.roi_range_margin;
roi_margin_mask = roi_margin.roi_margin_mask;
active_range_bins = roi_margin.active_range_bins;
roi_processing_count = roi_margin.roi_processing_count;
roi_extra_margin_bins = roi_margin.roi_extra_margin_bins;
print_roi_cfar_margin_selection(roi_margin, half_N, Nrx);

[sparse_range_doppler_rx, sparse_RD_power, sparse_fft_count, sparse_fft_work] = selective_doppler_processing_process(range_data_for_doppler, active_range_bins, half_N, Nchirp, Nrx, apply_clutter_removal);
print_selective_doppler_processing(active_range_bins, half_N, sparse_fft_count, sparse_fft_work);

radar_roi = struct("half_N",half_N,"Nchirp",Nchirp,"Tr",Tr,"Td",Td,"Gr",Gr,"Gd",Gd,"Pfa",Pfa,"Ntraining",Ntraining);
[roi_cfar_local, roi_noise_local, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, roi_stats, roi_det_range_idx, roi_det_doppler_idx, range_lookup] = roi_only_cfar_process(roi_bins, active_range_bins, sparse_RD_power, radar_roi);
print_roi_only_cfar(roi_stats);

regression = roi_cfar_regression(range_doppler_rx, active_range_bins, sparse_range_doppler_rx, RD_power, sparse_RD_power, cfar_map, roi_bins, roi_cfar_local, noise_map, roi_noise_local, peak_range_idx, peak_doppler_idx, roi_peak_range_idx, roi_peak_doppler_idx);
fft_difference = regression.fft_difference;
power_difference = regression.power_difference;
noise_difference = regression.noise_difference;
fft_match = regression.fft_match;
power_match = regression.power_match;
cfar_match = regression.cfar_match;
noise_match = regression.noise_match;
peak_match = regression.peak_match;



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


end
