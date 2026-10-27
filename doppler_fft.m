%% Doppler axis

PRF = 1 / Tchirp;

fd_axis = (-Nchirp/2:Nchirp/2-1) * PRF / Nchirp;

velocity_axis = fd_axis * lambda / 2;


%% Doppler Hann window

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);


%% Apply Doppler window

range_data_windowed = range_data_for_doppler .* reshape(doppler_window, 1, 1, Nchirp);


%% Doppler FFT for all RX channels

range_doppler_rx = fft(range_data_windowed, Nchirp, 3);

range_doppler_rx = fftshift(range_doppler_rx, 3);


%% Optional static clutter notch

if apply_clutter_removal

    center_bin = Nchirp / 2 + 1;

    notch_half_width = 1;

    range_doppler_rx(:, :, center_bin-notch_half_width:center_bin+notch_half_width) = 0;

end


%% Combine RX power for detection

RD_power = abs(range_doppler_rx(1:half_N, :, :)).^2;

RD_power = sum(RD_power, 2);

RD_power = squeeze(RD_power);

RD_mag = sqrt(RD_power);


%% Target range bins

[~, range_idx_1] = min(abs(range_axis - target_range(1)));

[~, range_idx_2] = min(abs(range_axis - target_range(2)));


%% Target 1 Doppler

doppler_target_1 = RD_mag(range_idx_1, :);

figure;

plot(velocity_axis, doppler_target_1);

xlabel("Radial Velocity [m/s]");

ylabel("Magnitude");

title("Target 1 Doppler FFT");

grid on;


%% Target 2 Doppler

doppler_target_2 = RD_mag(range_idx_2, :);

figure;

plot(velocity_axis, doppler_target_2);

xlabel("Radial Velocity [m/s]");

ylabel("Magnitude");

title("Target 2 Doppler FFT");

grid on;


%% Range-Doppler map

RD_db = 20*log10(RD_mag / max(RD_mag(:)) + eps);

figure;

imagesc(velocity_axis, range_axis, RD_db);

axis xy;

xlabel("Radial Velocity [m/s]");

ylabel("Range [m]");

title("Range-Doppler Map");

colorbar;

ylim([0 70]);

caxis([-60 0]);