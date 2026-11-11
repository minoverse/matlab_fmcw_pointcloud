
function [range_doppler_rx, RD_power, velocity_axis] = doppler_fft_process(range_data_for_doppler, cfg, apply_clutter_removal)

%% Doppler axis

PRF = 1 / cfg.Tchirp;

fd_axis = (-cfg.Nchirp/2:cfg.Nchirp/2-1) * PRF / cfg.Nchirp;

velocity_axis = fd_axis * cfg.lambda / 2;


%% Doppler Hann window

n_doppler = (0:cfg.Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/cfg.Nchirp);


%% Apply Doppler window

range_data_windowed = range_data_for_doppler .* ...
    reshape(doppler_window, 1, 1, cfg.Nchirp);


%% Doppler FFT

range_doppler_rx = fft(range_data_windowed, cfg.Nchirp, 3);

range_doppler_rx = fftshift(range_doppler_rx, 3);


%% Optional static clutter notch

if apply_clutter_removal

    center_bin = cfg.Nchirp / 2 + 1;

    notch_half_width = 1;

    range_doppler_rx(:, :, ...
        center_bin-notch_half_width:center_bin+notch_half_width) = 0;

end


%% Combine RX power

half_N = floor(cfg.Ns / 2);

RD_power = abs(range_doppler_rx(1:half_N, :, :)).^2;

RD_power = sum(RD_power, 2);

RD_power = squeeze(RD_power);

end
