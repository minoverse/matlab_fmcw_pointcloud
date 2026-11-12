function [sparse_range_doppler_rx, sparse_RD_power, sparse_fft_count, sparse_fft_work] = selective_doppler_processing_process(range_data_for_doppler, active_range_bins, half_N, Nchirp, Nrx, apply_clutter_removal)

%% ROI Bin 검증

active_range_bins = active_range_bins(:);

assert(~isempty(active_range_bins), "No active range bins");
assert(all(active_range_bins >= 1));
assert(all(active_range_bins <= half_N));
assert(all(active_range_bins == floor(active_range_bins)));
assert(length(unique(active_range_bins)) == length(active_range_bins));


%% 선택된 Range 데이터 추출

selected_range_data = range_data_for_doppler(active_range_bins, :, :);


%% Doppler Hann 윈도우 적용

n_doppler_sparse = (0:Nchirp-1)';

doppler_window_sparse = 0.5 - 0.5*cos(2*pi*n_doppler_sparse/Nchirp);

selected_windowed = selected_range_data .* reshape(doppler_window_sparse, 1, 1, Nchirp);


%% 선택적 Doppler FFT 수행

sparse_range_doppler_rx = fft(selected_windowed, Nchirp, 3);

sparse_range_doppler_rx = fftshift(sparse_range_doppler_rx, 3);


%% 선택적 정적 클러터 제거

if apply_clutter_removal

    center_bin = Nchirp / 2 + 1;

    notch_half_width = 1;

    sparse_range_doppler_rx(:, :, center_bin-notch_half_width:center_bin+notch_half_width) = 0;

end


%% RX 채널별 Power 합산

sparse_power_cube = abs(sparse_range_doppler_rx).^2;

sparse_power_sum = sum(sparse_power_cube, 2);

sparse_RD_power = reshape(sparse_power_sum, length(active_range_bins), Nchirp);


%% Doppler FFT 연산량 계산

sparse_fft_count = length(active_range_bins) * Nrx;

sparse_fft_work = sparse_fft_count * Nchirp * log2(Nchirp);


%% 결과 검증

assert(sparse_fft_count == length(active_range_bins) * Nrx);

end