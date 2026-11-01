
%% Phase 2 Group 4A Main

rng(1);

run("phase2_main.m");

addpath("tests");


%% ROI Merge

[roi_blocks, expanded_bins] = roi_merge_blocks(roi_bins, half_N, Tr+Gr);

assert(~isempty(roi_blocks));


%% Block Doppler Processing

[block_doppler_rx, block_RD_power, block_fft_count] = block_doppler_processing(range_data_for_doppler, roi_blocks, apply_clutter_removal);


%% ROI-Only CFAR

active_range_bins = expanded_bins;

sparse_RD_power = block_RD_power;

run("roi_only_cfar.m");


%% Regression

run("roi_block_regression.m");


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

