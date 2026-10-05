%% Range FFT for all chirps
% adc_cube dimensions:
% dim 1 = fast time samples
% dim 2 = RX channels
% dim 3 = chirps / slow time

range_data = fft(adc_cube, [], 1);


%% Doppler axis

PRF = 1 / Tchirp;

fd_axis = (-Nchirp/2 : Nchirp/2-1) * PRF / Nchirp;

velocity_axis = fd_axis * lambda / 2;


%% Range-Doppler FFT
% Doppler FFT along chirp dimension = dimension 3

range_doppler = fftshift( ...
    fft(range_data(:, 1, :), Nchirp, 3), ...
    3);

% Remove singleton RX dimension
range_doppler = squeeze(range_doppler);


%% Range axis

freq_axis = (0:Ns-1) * fs / Ns;

range_axis_full = c * freq_axis / (2 * S);


%% Keep positive range bins

half_N = floor(Ns / 2);

range_axis_rd = range_axis_full(1:half_N);

RD_mag = abs(range_doppler(1:half_N, :));


%% Find strongest Range-Doppler peak

[~, max_idx] = max(RD_mag(:));

[range_idx, doppler_idx] = ind2sub(size(RD_mag), max_idx);

estimated_range_rd = range_axis_rd(range_idx);
estimated_velocity = velocity_axis(doppler_idx);


%% Print result

fprintf("\n--- Range-Doppler Result ---\n");

fprintf("Ground Truth Range    = %.3f m\n", target_range);
fprintf("Estimated Range       = %.3f m\n", estimated_range_rd);

fprintf("Ground Truth Velocity = %.3f m/s\n", target_velocity);
fprintf("Estimated Velocity    = %.3f m/s\n", estimated_velocity);


%% Plot Doppler FFT at detected range bin

slow_signal = range_doppler(range_idx, :);

figure;
drawnow;
plot(velocity_axis, abs(slow_signal));

xlabel("Radial Velocity [m/s]");
ylabel("Magnitude");
title("Doppler FFT");

grid on;


%% Plot Range-Doppler Map

figure;
drawnow;
imagesc( ...
    velocity_axis, ...
    range_axis_rd, ...
    20 * log10(RD_mag + eps));

axis xy;

xlabel("Radial Velocity [m/s]");
ylabel("Range [m]");
title("Range-Doppler Map");

colorbar;

ylim([0 40]);