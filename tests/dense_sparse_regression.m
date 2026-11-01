
%% Dense vs Sparse Regression

%% Save Dense Reference

dense_doppler_reference = range_doppler_rx(1:half_N, :, :);

dense_power_reference = RD_power;

dense_cfar_reference = cfar_map;

dense_noise_reference = noise_map;

dense_det_range_reference = det_range_idx;
dense_det_doppler_reference = det_doppler_idx;

dense_peak_range_reference = peak_range_idx;
dense_peak_doppler_reference = peak_doppler_idx;
dense_peak_power_reference = peak_power;

dense_detection_count = nnz(cfar_map);


%% Full ROI

active_range_bins = (1:half_N)';

run("selective_doppler_processing.m");

assert(sparse_fft_count == half_N * Nrx);


%% Compare Complex Doppler FFT

fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_doppler_reference(:)));

fft_tolerance = 1e-10 * max(1, max(abs(dense_doppler_reference(:))));

fft_match = fft_difference <= fft_tolerance;

assert(fft_match, "Full ROI Doppler FFT mismatch");


%% Compare RD Power

power_difference = max(abs(sparse_RD_power(:) - dense_power_reference(:)));

power_tolerance = 1e-10 * max(1, max(abs(dense_power_reference(:))));

power_match = power_difference <= power_tolerance;

assert(power_match, "Full ROI RD power mismatch");


%% Run Original CFAR on Full ROI

RD_power = sparse_RD_power;

run("cfar2d.m");

run("peak_grouping.m");


%% Compare CFAR Detections

cfar_match = isequal(cfar_map, dense_cfar_reference);

assert(cfar_match, "Full ROI CFAR detection mismatch");


%% Compare Grouped Peaks

dense_peak_pairs = [dense_peak_range_reference(:), dense_peak_doppler_reference(:)];

sparse_peak_pairs = [peak_range_idx(:), peak_doppler_idx(:)];

peak_match = isequal(sortrows(dense_peak_pairs), sortrows(sparse_peak_pairs));

assert(peak_match, "Full ROI peak grouping mismatch");


%% Restore Dense Reference

RD_power = dense_power_reference;

cfar_map = dense_cfar_reference;

noise_map = dense_noise_reference;

det_range_idx = dense_det_range_reference;
det_doppler_idx = dense_det_doppler_reference;

peak_range_idx = dense_peak_range_reference;
peak_doppler_idx = dense_peak_doppler_reference;
peak_power = dense_peak_power_reference;


%% Results

fprintf("\n--- Dense vs Sparse Full ROI Regression ---\n");

fprintf("Full ROI bins = %d / %d\n", length(active_range_bins), half_N);

fprintf("Maximum FFT difference = %.3e\n", fft_difference);

fprintf("Maximum RD power difference = %.3e\n", power_difference);

fprintf("Dense CFAR detected cells = %d\n", dense_detection_count);

fprintf("Full ROI CFAR detected cells = %d\n", nnz(cfar_map));

fprintf("Complex Doppler FFT match = %d\n", fft_match);

fprintf("RD power match = %d\n", power_match);

fprintf("CFAR detection match = %d\n", cfar_match);

fprintf("Peak grouping match = %d\n", peak_match);

fprintf("Full ROI regression = PASS\n");
