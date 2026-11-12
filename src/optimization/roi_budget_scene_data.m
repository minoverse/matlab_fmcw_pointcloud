function data = roi_budget_scene_data(target_ranges, target_velocities, target_amplitudes, radar)

    %% 실험 Scene 설정

    target_ranges = target_ranges(:);
    target_velocities = target_velocities(:);
    target_amplitudes = target_amplitudes(:);

    num_targets = length(target_ranges);

    assert(length(target_velocities) == num_targets);
    assert(length(target_amplitudes) == num_targets);
    assert(num_targets > 0);


    %% ADC 신호 생성

    adc_cube = complex(zeros(radar.Ns, radar.Nrx, radar.Nchirp));

    t_fast = (0:radar.Ns-1)' / radar.fs;
    chirp_index = 0:radar.Nchirp-1;

    for k = 1:num_targets

        beat_frequency = radar.S * 2 * target_ranges(k) / radar.c;

        doppler_frequency = 2 * target_velocities(k) / radar.lambda;

        fast_wave = exp(1j*2*pi*beat_frequency*t_fast);

        slow_wave = exp(1j*2*pi*doppler_frequency*chirp_index*radar.Tchirp);

        base_signal = fast_wave * slow_wave;


        %% RX 채널별 공간 위상 계산

        if k == 1
            azimuth = 15;
            elevation = 10;
        elseif k == 2
            azimuth = -20;
            elevation = -8;
        else
            azimuth = -35 + 5*k;
            elevation = -10 + k;
        end

        direction_y = cosd(elevation)*sind(azimuth);
        direction_z = sind(elevation);

        for rx = 1:radar.Nrx

            spatial_phase = 2*pi/radar.lambda * (radar.rx_y(rx)*direction_y + radar.rx_z(rx)*direction_z);

            signal = target_amplitudes(k) * exp(1j*spatial_phase) * base_signal;

            adc_cube(:,rx,:) = adc_cube(:,rx,:) + reshape(signal,radar.Ns,1,radar.Nchirp);

        end

    end


    %% 복소 Gaussian Noise 추가

    noise_power = max(target_amplitudes)^2 / 10^(radar.snr_db/10);

    noise = sqrt(noise_power/2) * (randn(size(adc_cube)) + 1j*randn(size(adc_cube)));

    adc_cube = adc_cube + noise;


    %% Range FFT 수행

    n_range = (0:radar.Ns-1)';

    range_window = 0.5 - 0.5*cos(2*pi*n_range/(radar.Ns-1));

    adc_windowed = adc_cube .* reshape(range_window,radar.Ns,1,1);

    range_data = fft(adc_windowed,[],1);


    %% 비동기 합산 방식의 Range Energy 계산

    range_positive = range_data(1:radar.half_N,:,:);

    range_energy = squeeze(sum(sum(abs(range_positive).^2,3),2));


    %% 전체 Positive Range Bin에 Doppler FFT 적용

    n_doppler = (0:radar.Nchirp-1)';

    doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/radar.Nchirp);

    doppler_input = range_positive .* reshape(doppler_window,1,1,radar.Nchirp);

    dense_fft = fft(doppler_input,radar.Nchirp,3);

    dense_fft = fftshift(dense_fft,3);

    RD_power = reshape(sum(abs(dense_fft).^2,2),radar.half_N,radar.Nchirp);


    %% Dense CA-CFAR 파라미터

    Tr = radar.Tr;
    Td = radar.Td;
    Gr = radar.Gr;
    Gd = radar.Gd;

    range_margin = Tr+Gr;
    doppler_margin = Td+Gd;

    cfar_kernel = ones(2*range_margin+1,2*doppler_margin+1);

    cfar_kernel(range_margin+1-Gr:range_margin+1+Gr,doppler_margin+1-Gd:doppler_margin+1+Gd) = 0;

    Ntraining = nnz(cfar_kernel);

    assert(Ntraining == radar.Ntraining);

    alpha = Ntraining * (radar.Pfa^(-1/Ntraining)-1);


    %% Dense CA-CFAR Noise 및 Detection Map 계산

    training_sum = conv2(RD_power,cfar_kernel,"same");

    noise_map = zeros(radar.half_N,radar.Nchirp);

    cfar_map = false(radar.half_N,radar.Nchirp);

    valid_range = range_margin+1:radar.half_N-range_margin;

    valid_doppler = doppler_margin+1:radar.Nchirp-doppler_margin;

    noise_map(valid_range,valid_doppler) = training_sum(valid_range,valid_doppler)/Ntraining;

    cfar_map(valid_range,valid_doppler) = RD_power(valid_range,valid_doppler) > alpha*noise_map(valid_range,valid_doppler);


    %% Ground Truth 비교용 Range 및 Velocity 축 생성

    range_axis = radar.c*(0:radar.half_N-1)*radar.fs/(2*radar.S*radar.Ns);

    PRF = 1/radar.Tchirp;

    velocity_axis = (-radar.Nchirp/2:radar.Nchirp/2-1)*PRF*radar.lambda/(2*radar.Nchirp);


    %% 결과 구조체 반환

    data.range_data = range_data;
    data.dense_fft = dense_fft;
    data.RD_power = RD_power;
    data.cfar_map = cfar_map;
    data.noise_map = noise_map;
    data.range_energy = range_energy;

    data.range_axis = range_axis;
    data.velocity_axis = velocity_axis;

    data.target_ranges = target_ranges;
    data.target_velocities = target_velocities;

end