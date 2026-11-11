
function plot_doppler_fft(RD_power, range_axis, velocity_axis, cfg)

%% RD magnitude

RD_mag = sqrt(RD_power);


%% Target range bins

[~, range_idx_1] = min(abs(range_axis - cfg.target_range(1)));
[~, range_idx_2] = min(abs(range_axis - cfg.target_range(2)));


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

end
