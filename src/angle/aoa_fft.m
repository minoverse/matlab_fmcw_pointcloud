%% AoA estimation

Npeaks = length(peak_range_idx);

estimated_angles = zeros(Npeaks, 1);

estimated_snr = zeros(Npeaks, 1);

u_axis = (-NangleFFT/2:NangleFFT/2-1) * lambda / (NangleFFT * d_rx);

valid_angle = abs(u_axis) <= 1;

angle_axis = asind(u_axis(valid_angle));

for i = 1:Npeaks

    r = peak_range_idx(i);

    d = peak_doppler_idx(i);

    rx_vector = squeeze(range_doppler_rx(r, :, d));

    angle_fft = fft(rx_vector, NangleFFT);

    angle_fft = fftshift(angle_fft);

    angle_mag = abs(angle_fft(valid_angle));

    [~, angle_idx] = max(angle_mag);

    estimated_angles(i) = angle_axis(angle_idx);

    local_noise = noise_map(r, d);

    if local_noise > 0

        estimated_snr(i) = 10*log10(RD_power(r, d) / local_noise);

    else

        estimated_snr(i) = NaN;

    end

end

fprintf("\n--- AoA ---\n");

for i = 1:Npeaks

    fprintf("Detection %d: Angle = %.2f deg, SNR = %.2f dB\n", i, estimated_angles(i), estimated_snr(i));

end

%% Plot
figure;

stem(estimated_angles, ones(size(estimated_angles)));

xlabel("Angle [deg]");

ylabel("Detection");

title("Estimated Target Angles");

grid on;

xlim([-90 90]);

%% Save
if ~exist("plots/05_aoa", "dir")
    mkdir("plots/05_aoa");
end

exportgraphics(gcf, "plots/05_aoa/aoa_detections.png", "Resolution", 300);