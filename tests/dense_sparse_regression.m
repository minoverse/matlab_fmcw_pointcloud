function result = dense_sparse_regression(range_doppler_rx, RD_power, cfar_map, noise_map, det_range_idx, det_doppler_idx, peak_range_idx, peak_doppler_idx, peak_power, range_data_for_doppler, half_N, Nchirp, Nrx, apply_clutter_removal, cfg)
%% Dense vs Sparse Regression

%% 기존 Dense 결과 보존

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


%% 전체 Positive Range Bin 선택

active_range_bins = (1:half_N)';

[sparse_range_doppler_rx, sparse_RD_power, sparse_fft_count, sparse_fft_work] = selective_doppler_processing_process(range_data_for_doppler, active_range_bins, half_N, Nchirp, Nrx, apply_clutter_removal);

assert(sparse_fft_count == half_N * Nrx);


%% Complex Doppler FFT 비교

fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_doppler_reference(:)));

fft_tolerance = 1e-10 * max(1, max(abs(dense_doppler_reference(:))));

fft_match = fft_difference <= fft_tolerance;

assert(fft_match, "Full ROI Doppler FFT mismatch");


%% RD Power 비교

power_difference = max(abs(sparse_RD_power(:) - dense_power_reference(:)));

power_tolerance = 1e-10 * max(1, max(abs(dense_power_reference(:))));

power_match = power_difference <= power_tolerance;

assert(power_match, "Full ROI RD power mismatch");


%% 새 CFAR 함수로 Sparse 결과 검증

[sparse_cfar_map, sparse_noise_map, sparse_det_r, sparse_det_d, sparse_cfar_stats] = cfar2d_process(sparse_RD_power, cfg.cfar);


%% CFAR Detection 비교

cfar_match = isequal(sparse_cfar_map, dense_cfar_reference);

assert(cfar_match, "Full ROI CFAR detection mismatch");


%% 새 Peak Grouping 함수로 검증

[sparse_peak_r, sparse_peak_d, sparse_peak_power] = peak_grouping_process(sparse_RD_power, sparse_det_r, sparse_det_d);

dense_peak_pairs = [dense_peak_range_reference(:), dense_peak_doppler_reference(:)];
sparse_peak_pairs = [sparse_peak_r(:), sparse_peak_d(:)];

peak_match = isequal(sortrows(dense_peak_pairs), sortrows(sparse_peak_pairs));

assert(peak_match, "Full ROI peak grouping mismatch");


%% 기존 Dense 결과 유지

RD_power = dense_power_reference;
cfar_map = dense_cfar_reference;
noise_map = dense_noise_reference;

det_range_idx = dense_det_range_reference;
det_doppler_idx = dense_det_doppler_reference;

peak_range_idx = dense_peak_range_reference;
peak_doppler_idx = dense_peak_doppler_reference;
peak_power = dense_peak_power_reference;


%% 검증 결과 출력

fprintf("\n--- Dense vs Sparse Full ROI Regression ---\n");

fprintf("Full ROI bins = %d / %d\n", length(active_range_bins), half_N);
fprintf("Maximum FFT difference = %.3e\n", fft_difference);
fprintf("Maximum RD power difference = %.3e\n", power_difference);

fprintf("Dense CFAR detected cells = %d\n", dense_detection_count);
fprintf("Full ROI CFAR detected cells = %d\n", nnz(sparse_cfar_map));

fprintf("Complex Doppler FFT match = %d\n", fft_match);
fprintf("RD power match = %d\n", power_match);
fprintf("CFAR detection match = %d\n", cfar_match);
fprintf("Peak grouping match = %d\n", peak_match);

fprintf("Full ROI regression = PASS\n");

result = struct();
result.fft_difference = fft_difference;
result.power_difference = power_difference;
result.fft_match = fft_match;
result.power_match = power_match;
result.cfar_match = cfar_match;
result.peak_match = peak_match;
result.dense_detection_count = dense_detection_count;
end
