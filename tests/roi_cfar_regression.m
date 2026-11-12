function result = roi_cfar_regression(range_doppler_rx, active_range_bins, sparse_range_doppler_rx, RD_power, sparse_RD_power, cfar_map, roi_bins, roi_cfar_local, noise_map, roi_noise_local, peak_range_idx, peak_doppler_idx, roi_peak_range_idx, roi_peak_doppler_idx)

%% ROI CFAR Regression

%% Compare Complex Doppler FFT

dense_fft_selected = range_doppler_rx(active_range_bins, :, :);

fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_fft_selected(:)));

fft_tolerance = 1e-10 * max(1, max(abs(dense_fft_selected(:))));

fft_match = fft_difference <= fft_tolerance;

assert(fft_match, "Sparse Doppler FFT mismatch");


%% Compare RD Power

dense_power_selected = RD_power(active_range_bins, :);

power_difference = max(abs(sparse_RD_power(:) - dense_power_selected(:)));

power_tolerance = 1e-10 * max(1, max(abs(dense_power_selected(:))));

power_match = power_difference <= power_tolerance;

assert(power_match, "Sparse RD power mismatch");


%% Compare CFAR Only at ROI CUTs

dense_cfar_roi = cfar_map(roi_bins, :);

cfar_match = isequal(roi_cfar_local, dense_cfar_roi);

assert(cfar_match, "ROI CFAR mismatch");


%% Compare CFAR Noise Estimates

dense_noise_roi = noise_map(roi_bins, :);

noise_difference = max(abs(roi_noise_local(:) - dense_noise_roi(:)));

noise_tolerance = 1e-10 * max(1, max(abs(dense_noise_roi(:))));

noise_match = noise_difference <= noise_tolerance;

assert(noise_match, "CFAR noise estimate mismatch");


%% Compare Grouped Peaks Within ROI

dense_peak_pairs = [peak_range_idx(:), peak_doppler_idx(:)];

dense_peak_pairs = dense_peak_pairs(ismember(dense_peak_pairs(:, 1), roi_bins), :);

roi_peak_pairs = [roi_peak_range_idx(:), roi_peak_doppler_idx(:)];

peak_match = isequal(sortrows(dense_peak_pairs), sortrows(roi_peak_pairs));

assert(peak_match, "ROI peak grouping mismatch");


%% Results

fprintf("\n--- ROI CFAR Regression ---\n");

fprintf("ROI bins tested = %d\n", length(roi_bins));

fprintf("Processed bins including margin = %d\n", length(active_range_bins));

fprintf("Maximum FFT difference = %.3e\n", fft_difference);

fprintf("Maximum RD power difference = %.3e\n", power_difference);

fprintf("Maximum CFAR noise difference = %.3e\n", noise_difference);

fprintf("FFT match = %d\n", fft_match);

fprintf("RD power match = %d\n", power_match);

fprintf("CFAR detection match = %d\n", cfar_match);

fprintf("CFAR noise match = %d\n", noise_match);

fprintf("Peak grouping match = %d\n", peak_match);

fprintf("ROI CFAR regression = PASS\n");


result = struct();
result.fft_difference = fft_difference;
result.power_difference = power_difference;
result.noise_difference = noise_difference;
result.fft_match = fft_match;
result.power_match = power_match;
result.cfar_match = cfar_match;
result.noise_match = noise_match;
result.peak_match = peak_match;
end
