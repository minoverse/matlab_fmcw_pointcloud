%% Azimuth and Elevation estimation

Npeaks = length(peak_range_idx);

estimated_azimuth = zeros(Npeaks, 1);

estimated_elevation = zeros(Npeaks, 1);

estimated_snr = zeros(Npeaks, 1);


%% Spatial frequency axes

u_az_axis = (-NangleFFT_az/2:NangleFFT_az/2-1) * lambda / (NangleFFT_az * d_az);

u_el_axis = (-NangleFFT_el/2:NangleFFT_el/2-1) * lambda / (NangleFFT_el * d_el);


[u_az_grid, u_el_grid] = meshgrid(u_az_axis, u_el_axis);

physical_mask = u_az_grid.^2 + u_el_grid.^2 <= 1;


%% Process each detected target

for i = 1:Npeaks

    r = peak_range_idx(i);

    d = peak_doppler_idx(i);


    rx_vector = squeeze(range_doppler_rx(r, :, d));


    %% Map RX channels into 2 x 2 planar array

    rx_matrix = zeros(2, 2);

    rx_matrix(1, 1) = rx_vector(1);

    rx_matrix(1, 2) = rx_vector(2);

    rx_matrix(2, 1) = rx_vector(3);

    rx_matrix(2, 2) = rx_vector(4);


    %% 2D angle FFT

    angle_fft = fft2(rx_matrix, NangleFFT_el, NangleFFT_az);

    angle_fft = fftshift(angle_fft);

    angle_power = abs(angle_fft).^2;


    %% Remove non-physical direction cosines

    angle_power(~physical_mask) = 0;


    %% Find strongest angle peak

    [~, peak_index] = max(angle_power(:));

    [el_idx, az_idx] = ind2sub(size(angle_power), peak_index);


    direction_y = u_az_axis(az_idx);

    direction_z = u_el_axis(el_idx);


    %% Elevation

    elevation = asind(direction_z);


    %% Azimuth

    cos_elevation = cosd(elevation);

    azimuth_argument = direction_y / cos_elevation;

    azimuth_argument = max(-1, min(1, azimuth_argument));

    azimuth = asind(azimuth_argument);


    estimated_azimuth(i) = azimuth;

    estimated_elevation(i) = elevation;


    %% SNR estimate

    local_noise = noise_map(r, d);

    if local_noise > 0

        estimated_snr(i) = 10*log10(RD_power(r, d) / local_noise);

    else

        estimated_snr(i) = NaN;

    end

end


%% Print

fprintf("\n--- Azimuth / Elevation AoA ---\n");

for i = 1:Npeaks

    fprintf("Detection %d: Azimuth = %.2f deg, Elevation = %.2f deg, SNR = %.2f dB\n", i, estimated_azimuth(i), estimated_elevation(i), estimated_snr(i));

end


%% Plot

figure;

scatter(estimated_azimuth, estimated_elevation, 80, "filled");

xlabel("Azimuth [deg]");

ylabel("Elevation [deg]");

title("Azimuth / Elevation Detections");

grid on;

xlim([-90 90]);

ylim([-60 60]);


%% Save

if ~exist("plots/05_aoa", "dir")

    mkdir("plots/05_aoa");

end


exportgraphics(gcf, "plots/05_aoa/aoa_detections.png", "Resolution", 300);