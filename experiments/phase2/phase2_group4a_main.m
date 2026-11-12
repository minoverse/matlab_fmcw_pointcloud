function phase2_group4a_main()

%% Phase 2 Group 4A Main

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


%% ROI Merge

[roi_blocks, expanded_bins] = roi_merge_blocks(roi_bins, half_N, Tr+Gr);

assert(~isempty(roi_blocks));


%% Block Doppler Processing

[block_doppler_rx, block_RD_power, block_fft_count] = block_doppler_processing(range_data_for_doppler, roi_blocks, apply_clutter_removal);


%% ROI-Only CFAR

active_range_bins = expanded_bins;

sparse_RD_power = block_RD_power;

radar_roi = struct("half_N",half_N,"Nchirp",Nchirp,"Tr",Tr,"Td",Td,"Gr",Gr,"Gd",Gd,"Pfa",Pfa,"Ntraining",Ntraining);
[roi_cfar_local, roi_noise_local, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, roi_stats, roi_det_range_idx, roi_det_doppler_idx, range_lookup] = roi_only_cfar_process(roi_bins, active_range_bins, sparse_RD_power, radar_roi);
print_roi_only_cfar(roi_stats);


%% Regression

regression = roi_block_regression(roi_bins, active_range_bins, roi_cfar_local, roi_noise_local, roi_det_range_idx, roi_det_doppler_idx, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, range_lookup, half_N, Tr, Gr, expanded_bins, roi_blocks, range_doppler_rx, block_doppler_rx, RD_power, block_RD_power, cfar_map, noise_map, peak_range_idx, peak_doppler_idx, block_fft_count, Nrx, Nchirp, Td, Gd, range_data_for_doppler, apply_clutter_removal, Pfa, Ntraining);
selected_cfar_cuts = regression.selected_cfar_cuts;
selected_cfar_map = regression.selected_cfar_map;
selected_noise_map = regression.selected_noise_map;
selected_det_range_idx = regression.selected_det_range_idx;
selected_det_doppler_idx = regression.selected_det_doppler_idx;
selected_fft_error = regression.selected_fft_error;
selected_power_error = regression.selected_power_error;
selected_noise_error = regression.selected_noise_error;
full_fft_error = regression.full_fft_error;
full_power_error = regression.full_power_error;
selected_roi_bins = regression.selected_roi_bins;
selected_active_bins = regression.selected_active_bins;
selected_peak_range_idx = regression.selected_peak_range_idx;
selected_peak_doppler_idx = regression.selected_peak_doppler_idx;



%% Workload Counters

num_blocks = size(roi_blocks, 1);

merged_bin_count = length(expanded_bins);

unmerged_bin_slots = 0;

for i = 1:length(roi_bins)

    first_bin = max(1, roi_bins(i) - Tr - Gr);
    last_bin = min(half_N, roi_bins(i) + Tr + Gr);

    unmerged_bin_slots = unmerged_bin_slots + last_bin - first_bin + 1;

end

duplicate_slots_removed = unmerged_bin_slots - merged_bin_count;

fair_dense_fft_count = half_N * Nrx;

fft_reduction = 100 * (1 - block_fft_count / fair_dense_fft_count);

cfar_reduction = 100 * (1 - selected_cfar_cuts / cfar_cut_count);


%% Logical Buffer Size Proxies

complex_bytes = 16;

block_lengths = roi_blocks(:,2) - roi_blocks(:,1) + 1;

max_block_bins = max(block_lengths);

max_block_input_bytes = max_block_bins * Nrx * Nchirp * complex_bytes;

collected_fft_bytes = merged_bin_count * Nrx * Nchirp * complex_bytes;


%% Validation

assert(block_fft_count == merged_bin_count * Nrx);

assert(duplicate_slots_removed >= 0);

assert(merged_bin_count <= half_N);

assert(selected_cfar_cuts <= cfar_cut_count);


%% Output Folder

output_dir = "output/phase2/08_roi_merge_block_processing";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Save Block CSV

block_ids = (1:num_blocks)';

block_start = roi_blocks(:,1);

block_end = roi_blocks(:,2);

block_start_range = range_axis(block_start)';

block_end_range = range_axis(block_end)';

block_table = table(block_ids, block_start, block_end, block_lengths, block_start_range, block_end_range);

block_table.Properties.VariableNames = {'BlockID', 'StartBin', 'EndBin', 'BinCount', 'StartRange_m', 'EndRange_m'};

csv_file = fullfile(output_dir, "merged_roi_processing_blocks.csv");

writetable(block_table, csv_file);


%% Plot ROI Blocks

processing_mask = false(half_N, 1);

processing_mask(expanded_bins) = true;

fig = figure("Position", [100 100 1100 750]);

subplot(2,1,1);

stairs(range_axis, double(processing_mask), "b", "LineWidth", 1.4);

hold on;

stem(range_axis(roi_bins), 1.15 * ones(length(roi_bins),1), "r", "filled");

xlabel("Range [m]");
ylabel("Selected");

title("Merged ROI Blocks and Original CFAR CUT Bins");

legend("Merged Processing Blocks", "Original ROI Bins");

xlim([0 70]);
ylim([0 1.4]);
grid on;


%% Plot Detection Regression

subplot(2,1,2);

imagesc(velocity_axis, range_axis, RD_db);

axis xy;
hold on;

[dense_r, dense_d] = find(cfar_map);

h1 = plot(velocity_axis(dense_d), range_axis(dense_r), "ro", "MarkerSize", 8);

h2 = plot(velocity_axis(selected_det_doppler_idx), range_axis(selected_det_range_idx), "g+", "MarkerSize", 9, "LineWidth", 1.5);

xlabel("Velocity [m/s]");
ylabel("Range [m]");

title("Dense CFAR vs ROI Block CFAR");

legend([h1, h2], "Dense CFAR", "ROI Block CFAR");

ylim([0 70]);
caxis([-60 0]);
colorbar;


%% Save PNG

plot_file = fullfile(output_dir, "merged_roi_blocks_and_cfar_regression.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save TXT

txt_file = fullfile(output_dir, "roi_block_processing_validation.txt");

file = fopen(txt_file, "w");

assert(file ~= -1, "Cannot create Group 4A result file");

fprintf(file, "=== Phase 2 Group 4A: ROI Merge and Block Processing ===\n\n");

fprintf(file, "[ROI Merge]\n");
fprintf(file, "Original ROI bins = %d\n", length(roi_bins));
fprintf(file, "CFAR range margin = %d\n", Tr+Gr);
fprintf(file, "Merged block count = %d\n", num_blocks);
fprintf(file, "Naive expanded bin slots = %d\n", unmerged_bin_slots);
fprintf(file, "Unique processing bins = %d\n", merged_bin_count);
fprintf(file, "Duplicate slots removed vs naive = %d\n", duplicate_slots_removed);

fprintf(file, "\n[Block Doppler]\n");
fprintf(file, "Fair Dense FFT count = %d\n", fair_dense_fft_count);
fprintf(file, "Block Doppler FFT count = %d\n", block_fft_count);
fprintf(file, "Doppler FFT workload reduction = %.2f%%\n", fft_reduction);

fprintf(file, "\n[ROI CFAR]\n");
fprintf(file, "Dense CFAR CUTs = %d\n", cfar_cut_count);
fprintf(file, "ROI CFAR CUTs = %d\n", selected_cfar_cuts);
fprintf(file, "CFAR CUT workload reduction = %.2f%%\n", cfar_reduction);
fprintf(file, "Dense ROI detections = %d\n", nnz(cfar_map(roi_bins,:)));
fprintf(file, "Block ROI detections = %d\n", nnz(selected_cfar_map));

fprintf(file, "\n[Logical Buffer Proxies]\n");
fprintf(file, "Largest block bins = %d\n", max_block_bins);
fprintf(file, "Largest block input = %.4f MiB\n", max_block_input_bytes / 2^20);
fprintf(file, "Collected FFT output = %.4f MiB\n", collected_fft_bytes / 2^20);

fprintf(file, "\n[Regression]\n");
fprintf(file, "Selected FFT error = %.3e\n", selected_fft_error);
fprintf(file, "Selected RD power error = %.3e\n", selected_power_error);
fprintf(file, "Selected CFAR noise error = %.3e\n", selected_noise_error);
fprintf(file, "Full ROI FFT error = %.3e\n", full_fft_error);
fprintf(file, "Full ROI RD power error = %.3e\n", full_power_error);
fprintf(file, "Overlapping ROI merge = PASS\n");
fprintf(file, "Adjacent ROI merge = PASS\n");
fprintf(file, "Edge ROI merge = PASS\n");
fprintf(file, "Empty ROI handling = PASS\n");
fprintf(file, "Full ROI regression = PASS\n");
fprintf(file, "Selected ROI regression = PASS\n");
fprintf(file, "Overall = PASS\n");

fprintf(file, "\n[Limitations]\n");
fprintf(file, "Group 3A already removes duplicate bins using a mask.\n");
fprintf(file, "Block merging is not an additional FFT saving over Group 3A.\n");
fprintf(file, "FFT processing is grouped by contiguous range block.\n");
fprintf(file, "All block outputs are still collected for ROI CFAR.\n");
fprintf(file, "Buffer sizes are logical proxies, not measured peak RAM.\n");
fprintf(file, "No CPU latency, DMA or power improvement is claimed.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 4A Validation ---\n");

fprintf("Original ROI bins = %d\n", length(roi_bins));

fprintf("Merged ROI blocks = %d\n", num_blocks);

fprintf("Unique expanded range bins = %d\n", merged_bin_count);

fprintf("Duplicate slots removed vs naive = %d\n", duplicate_slots_removed);

fprintf("Block Doppler FFTs = %d\n", block_fft_count);

fprintf("Doppler FFT workload reduction = %.2f%%\n", fft_reduction);

fprintf("ROI CFAR CUTs = %d\n", selected_cfar_cuts);

fprintf("ROI Block regression = PASS\n");

fprintf("Full ROI regression = PASS\n");

fprintf("Group 4A validation passed.\n");

fprintf("CSV saved: %s\n", csv_file);
fprintf("PNG saved: %s\n", plot_file);
fprintf("TXT saved: %s\n", txt_file);


end
