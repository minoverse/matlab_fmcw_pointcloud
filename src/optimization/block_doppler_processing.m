
function [block_doppler_rx, block_RD_power, block_fft_count] = block_doppler_processing(range_data_for_doppler, blocks, apply_clutter_removal)

%% Dimensions

Nrx = size(range_data_for_doppler, 2);
Nchirp = size(range_data_for_doppler, 3);

total_bins = sum(blocks(:, 2) - blocks(:, 1) + 1);


%% Allocate Output

block_doppler_rx = complex(zeros(total_bins, Nrx, Nchirp));

block_RD_power = zeros(total_bins, Nchirp);


%% Doppler Hann Window

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);


%% Process Each Block

position = 1;

for i = 1:size(blocks, 1)

    first_bin = blocks(i, 1);
    last_bin = blocks(i, 2);

    count = last_bin - first_bin + 1;


    %% Select Contiguous Range Block

    block_input = range_data_for_doppler(first_bin:last_bin, :, :);

    block_windowed = block_input .* reshape(doppler_window, 1, 1, Nchirp);


    %% Doppler FFT

    fft_data = fft(block_windowed, Nchirp, 3);

    fft_data = fftshift(fft_data, 3);


    %% Optional Static Clutter Notch

    if apply_clutter_removal

        center_bin = floor(Nchirp / 2) + 1;

        notch_half_width = 1;

        fft_data(:, :, center_bin-notch_half_width:center_bin+notch_half_width) = 0;

    end


    %% Store Complex FFT and RX Power

    output_rows = position:position+count-1;

    block_doppler_rx(output_rows, :, :) = fft_data;

    block_RD_power(output_rows, :) = reshape(sum(abs(fft_data).^2, 2), count, Nchirp);

    position = position + count;

end


%% Workload Counter

block_fft_count = total_bins * Nrx;

assert(position == total_bins + 1);
assert(size(block_RD_power, 1) == total_bins);

end

