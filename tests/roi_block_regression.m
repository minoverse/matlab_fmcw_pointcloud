
%% Group 4A - ROI Block Regression

%% Save Selected ROI Results

selected_roi_bins = roi_bins;
selected_active_bins = active_range_bins;

selected_cfar_map = roi_cfar_local;
selected_noise_map = roi_noise_local;

selected_det_range_idx = roi_det_range_idx;
selected_det_doppler_idx = roi_det_doppler_idx;

selected_peak_range_idx = roi_peak_range_idx;
selected_peak_doppler_idx = roi_peak_doppler_idx;

selected_cfar_cuts = roi_cfar_cut_count;
selected_range_lookup = range_lookup;


%% Check Merged Range Coverage

expected_mask = false(half_N, 1);

for i = 1:length(selected_roi_bins)

    r = selected_roi_bins(i);

    first_bin = max(1, r - Tr - Gr);
    last_bin = min(half_N, r + Tr + Gr);

    expected_mask(first_bin:last_bin) = true;

end

assert(isequal(expanded_bins, find(expected_mask)));

assert(sum(roi_blocks(:,2) - roi_blocks(:,1) + 1) == length(expanded_bins));


%% Compare Block FFT Against Dense

dense_selected_fft = range_doppler_rx(expanded_bins, :, :);

selected_fft_error = max(abs(block_doppler_rx(:) - dense_selected_fft(:)));

selected_fft_tolerance = 1e-10 * max(1, max(abs(dense_selected_fft(:))));

assert(selected_fft_error <= selected_fft_tolerance, "Block Doppler FFT mismatch");


%% Compare Range-Doppler Power

dense_selected_power = RD_power(expanded_bins, :);

selected_power_error = max(abs(block_RD_power(:) - dense_selected_power(:)));

selected_power_tolerance = 1e-10 * max(1, max(abs(dense_selected_power(:))));

assert(selected_power_error <= selected_power_tolerance, "Block RD power mismatch");


%% Compare ROI-Only CFAR Against Dense ROI CUTs

dense_roi_cfar = cfar_map(selected_roi_bins, :);

assert(isequal(selected_cfar_map, dense_roi_cfar), "ROI CFAR mismatch");

dense_roi_noise = noise_map(selected_roi_bins, :);

selected_noise_error = max(abs(selected_noise_map(:) - dense_roi_noise(:)));

selected_noise_tolerance = 1e-10 * max(1, max(abs(dense_roi_noise(:))));

assert(selected_noise_error <= selected_noise_tolerance);


%% Compare Peak Grouping

dense_peak_pairs = [peak_range_idx(:), peak_doppler_idx(:)];

dense_peak_pairs = dense_peak_pairs(ismember(dense_peak_pairs(:,1), selected_roi_bins), :);

selected_peak_pairs = [selected_peak_range_idx(:), selected_peak_doppler_idx(:)];

assert(isequal(sortrows(dense_peak_pairs), sortrows(selected_peak_pairs)));


%% Check Workload Counters

expected_cfar_cuts = nnz(selected_roi_bins > Tr+Gr & selected_roi_bins <= half_N-Tr-Gr) * (Nchirp - 2*(Td+Gd));

assert(block_fft_count == length(expanded_bins) * Nrx);
assert(selected_cfar_cuts == expected_cfar_cuts);


%% Synthetic Merge Tests

[test_blocks, test_bins] = roi_merge_blocks([30, 32, 70], half_N, 10);

assert(isequal(test_blocks, [20, 42; 60, 80]));
assert(length(test_bins) == 44);


[test_blocks, test_bins] = roi_merge_blocks([30, 51], half_N, 10);

assert(isequal(test_blocks, [20, 61]));
assert(length(test_bins) == 42);


[test_blocks, test_bins] = roi_merge_blocks([1, half_N], half_N, 10);

assert(isequal(test_blocks, [1, 11; half_N-10, half_N]));
assert(length(test_bins) == 22);


[test_blocks, test_bins] = roi_merge_blocks([], half_N, 10);

assert(isempty(test_blocks));
assert(isempty(test_bins));


[test_blocks, test_bins] = roi_merge_blocks((1:half_N)', half_N, 10);

assert(isequal(test_blocks, [1, half_N]));
assert(length(test_bins) == half_N);


%% Empty ROI Block FFT Test

[empty_fft, empty_power, empty_count] = block_doppler_processing(range_data_for_doppler, zeros(0,2), false);

assert(size(empty_fft,1) == 0);
assert(size(empty_power,1) == 0);
assert(empty_count == 0);


%% Full ROI Regression

[full_blocks, full_bins] = roi_merge_blocks((1:half_N)', half_N, Tr+Gr);

[full_fft, full_power, full_fft_count] = block_doppler_processing(range_data_for_doppler, full_blocks, apply_clutter_removal);

assert(full_fft_count == half_N * Nrx);


%% Full ROI Complex FFT

dense_full_fft = range_doppler_rx(1:half_N, :, :);

full_fft_error = max(abs(full_fft(:) - dense_full_fft(:)));

full_fft_tolerance = 1e-10 * max(1, max(abs(dense_full_fft(:))));

assert(full_fft_error <= full_fft_tolerance);


%% Full ROI Power

full_power_error = max(abs(full_power(:) - RD_power(:)));

full_power_tolerance = 1e-10 * max(1, max(abs(RD_power(:))));

assert(full_power_error <= full_power_tolerance);


%% Full ROI CFAR and Peak Grouping

roi_bins = (1:half_N)';

active_range_bins = full_bins;

sparse_RD_power = full_power;

run("roi_only_cfar.m");

assert(isequal(roi_cfar_local, cfar_map));
assert(roi_cfar_cut_count == 25920);

full_noise_error = max(abs(roi_noise_local(:) - noise_map(:)));

full_noise_tolerance = 1e-10 * max(1, max(abs(noise_map(:))));

assert(full_noise_error <= full_noise_tolerance);

full_peak_pairs = [roi_peak_range_idx(:), roi_peak_doppler_idx(:)];

dense_full_peaks = [peak_range_idx(:), peak_doppler_idx(:)];

assert(isequal(sortrows(full_peak_pairs), sortrows(dense_full_peaks)));


%% Restore Selected ROI Results

roi_bins = selected_roi_bins;
active_range_bins = selected_active_bins;

sparse_RD_power = block_RD_power;

roi_cfar_local = selected_cfar_map;
roi_noise_local = selected_noise_map;

roi_det_range_idx = selected_det_range_idx;
roi_det_doppler_idx = selected_det_doppler_idx;

roi_peak_range_idx = selected_peak_range_idx;
roi_peak_doppler_idx = selected_peak_doppler_idx;

roi_cfar_cut_count = selected_cfar_cuts;
range_lookup = selected_range_lookup;


%% Results

fprintf("\n--- Group 4A ROI Block Regression ---\n");

fprintf("Selected ROI FFT error = %.3e\n", selected_fft_error);
fprintf("Selected ROI power error = %.3e\n", selected_power_error);
fprintf("Selected ROI CFAR noise error = %.3e\n", selected_noise_error);

fprintf("Synthetic merge tests = PASS\n");
fprintf("Selected ROI CFAR / peak regression = PASS\n");

fprintf("Full ROI FFT error = %.3e\n", full_fft_error);
fprintf("Full ROI power error = %.3e\n", full_power_error);
fprintf("Full ROI CFAR / peak regression = PASS\n");

fprintf("Group 4A regression = PASS\n");

