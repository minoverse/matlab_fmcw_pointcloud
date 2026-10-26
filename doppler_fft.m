%% Doppler axis

PRF = 1 / Tchirp;

fd_axis = ...
    (-Nchirp/2:Nchirp/2-1) * PRF / Nchirp;

velocity_axis = fd_axis * lambda / 2;


%% Doppler Hann window

n_doppler = (0:Nchirp-1)';

doppler_window = ...
    0.5 - 0.5*cos(2*pi*n_doppler/(Nchirp-1));


%% Apply Doppler window

range_data_windowed = ...
    range_data .* reshape(doppler_window, 1, 1, Nchirp);


%% Doppler FFT

range_doppler = fftshift( ...
    fft(range_data_windowed(:, 1, :), Nchirp, 3), ...
    3);

range_doppler = squeeze(range_doppler);


%% Range-Doppler magnitude

RD_mag = abs(range_doppler(1:half_N, :));


%% Target 1 range bin

[~, range_idx_1] = ...
    min(abs(range_axis - target_range(1)));

doppler_target_1 = ...
    RD_mag(range_idx_1, :);


%% Target 2 range bin

[~, range_idx_2] = ...
    min(abs(range_axis - target_range(2)));

doppler_target_2 = ...
    RD_mag(range_idx_2, :);


%% Plot Target 1 Doppler

figure;

plot(velocity_axis, doppler_target_1);

xlabel("Radial Velocity [m/s]");
ylabel("Magnitude");

title("Target 1 Doppler FFT - Hann Window");

grid on;


%% Plot Target 2 Doppler

figure;

plot(velocity_axis, doppler_target_2);

xlabel("Radial Velocity [m/s]");
ylabel("Magnitude");

title("Target 2 Doppler FFT - Hann Window");

grid on;


%% Range-Doppler Map

RD_db = 20*log10(RD_mag / max(RD_mag(:)) + eps);

figure;

imagesc( ...
    velocity_axis, ...
    range_axis, ...
    RD_db);

axis xy;

xlabel("Radial Velocity [m/s]");
ylabel("Range [m]");

title("Range-Doppler Map - Hann Window");

colorbar;

ylim([0 70]);

caxis([-60 0]);