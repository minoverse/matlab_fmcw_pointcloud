
function [estimated_azimuth, estimated_elevation, estimated_snr] = aoa_fft_process(range_doppler_rx, RD_power, noise_map, peak_range_idx, peak_doppler_idx, cfg)

    %% 결과 배열 초기화

    Npeaks = length(peak_range_idx);

    estimated_azimuth = zeros(Npeaks, 1);
    estimated_elevation = zeros(Npeaks, 1);
    estimated_snr = zeros(Npeaks, 1);


    %% 공간 주파수 축 계산

    u_az_axis = (-cfg.NangleFFT_az/2:cfg.NangleFFT_az/2-1) * cfg.lambda / (cfg.NangleFFT_az * cfg.d_az);

    u_el_axis = (-cfg.NangleFFT_el/2:cfg.NangleFFT_el/2-1) * cfg.lambda / (cfg.NangleFFT_el * cfg.d_el);

    [u_az_grid, u_el_grid] = meshgrid(u_az_axis, u_el_axis);

    % 물리적으로 가능한 방향만 선택
    physical_mask = u_az_grid.^2 + u_el_grid.^2 <= 1;


    %% 검출된 Peak별 AoA 추정

    for i = 1:Npeaks

        r = peak_range_idx(i);
        d = peak_doppler_idx(i);

        % 해당 Range Doppler Bin의 RX 신호 추출
        rx_vector = squeeze(range_doppler_rx(r, :, d));


        %% RX 4개 -> 2x2 평면 배열 배치

        rx_matrix = zeros(2, 2);

        rx_matrix(1, 1) = rx_vector(1);
        rx_matrix(1, 2) = rx_vector(2);
        rx_matrix(2, 1) = rx_vector(3);
        rx_matrix(2, 2) = rx_vector(4);


        %% 2D Angle FFT

        angle_fft = fft2(rx_matrix, cfg.NangleFFT_el, cfg.NangleFFT_az);

        angle_fft = fftshift(angle_fft);

        angle_power = abs(angle_fft).^2;


        %% 불가능한 방향 제거

        angle_power(~physical_mask) = 0;


        %% Angle FFT에서 최대 Power 위치 탐색

        [~, peak_index] = max(angle_power(:));

        [el_idx, az_idx] = ind2sub(size(angle_power), peak_index);

        direction_y = u_az_axis(az_idx);
        direction_z = u_el_axis(el_idx);


        %% Elevation 계산

        elevation = asind(direction_z);


        %% Azimuth 계산

        cos_elevation = cosd(elevation);

        azimuth_argument = direction_y / cos_elevation;

        azimuth_argument = max(-1, min(1, azimuth_argument));

        azimuth = asind(azimuth_argument);


        %% 추정 결과 저장

        estimated_azimuth(i) = azimuth;
        estimated_elevation(i) = elevation;


        %% CFAR Noise 기반 SNR 추정

        local_noise = noise_map(r, d);

        if local_noise > 0

            estimated_snr(i) = 10*log10(RD_power(r, d) / local_noise);

        else

            estimated_snr(i) = NaN;

        end

    end

end
