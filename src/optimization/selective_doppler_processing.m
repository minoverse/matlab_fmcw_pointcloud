
%% Selective Doppler Processing

%% Check ROI Bins

active_range_bins = active_range_bins(:);

assert(~isempty(active_range_bins), "No active range bins");
assert(all(active_range_bins >= 1));
assert(all(active_range_bins <= half_N));
assert(all(active_range_bins == floor(active_range_bins)));
assert(length(unique(active_range_bins)) == length(active_range_bins));


%% Select Range Data

selected_range_data = range_data_for_doppler(active_range_bins, :, :);


%% Doppler Hann Window

n_doppler_sparse = (0:Nchirp-1)';

doppler_window_sparse = 0.5 - 0.5*cos(2*pi*n_doppler_sparse/Nchirp);

selected_windowed = selected_range_data .* reshape(doppler_window_sparse, 1, 1, Nchirp);


%% Selective Doppler FFT

sparse_range_doppler_rx = fft(selected_windowed, Nchirp, 3);

sparse_range_doppler_rx = fftshift(sparse_range_doppler_rx, 3);


%% Optional Static Clutter Notch

if apply_clutter_removal

    center_bin = Nchirp / 2 + 1;

    notch_half_width = 1;

    sparse_range_doppler_rx(:, :, center_bin-notch_half_width:center_bin+notch_half_width) = 0;

end


%% Combine RX Power

sparse_power_cube = abs(sparse_range_doppler_rx).^2;

sparse_power_sum = sum(sparse_power_cube, 2);

sparse_RD_power = reshape(sparse_power_sum, length(active_range_bins), Nchirp);


%% Workload Counter

sparse_fft_count = length(active_range_bins) * Nrx;

sparse_fft_work = sparse_fft_count * Nchirp * log2(Nchirp);


%% Validation

assert(sparse_fft_count == length(active_range_bins) * Nrx);


%% Results

fprintf("\n--- Selective Doppler Processing ---\n");

fprintf("Active range bins = %d / %d\n", length(active_range_bins), half_N);

fprintf("Doppler FFTs executed = %d\n", sparse_fft_count);

fprintf("Doppler FFT work proxy = %.0f\n", sparse_fft_work);

