
function [adc_cube, t_fast, noise] = simulate_adc_cube_process(cfg, snr_db)

%% Fast-time axis

n = 0:cfg.Ns-1;

t_fast = n / cfg.fs;


%% ADC cube

adc_cube = complex(zeros(cfg.Ns, cfg.Nrx, cfg.Nchirp));


%% Target signal generation

for m = 0:cfg.Nchirp-1

    for rx = 1:cfg.Nrx

        rx_signal = complex(zeros(1, cfg.Ns));

        for k = 1:cfg.Ntargets

            fast_phase = 2*pi * cfg.fb(k) * t_fast;

            slow_phase = 2*pi * cfg.fd(k) * m * cfg.Tchirp;

            azimuth = cfg.target_angle(k);
            elevation = cfg.target_elevation(k);

            direction_y = cosd(elevation) * sind(azimuth);
            direction_z = sind(elevation);

            spatial_phase = 2*pi/cfg.lambda * ...
                (cfg.rx_y(rx) * direction_y + ...
                cfg.rx_z(rx) * direction_z);

            target_signal = cfg.target_amplitude(k) * ...
                exp(1j * (fast_phase + slow_phase + spatial_phase));

            rx_signal = rx_signal + target_signal;

        end

        adc_cube(:, rx, m+1) = rx_signal.';

    end

end


%% Complex Gaussian noise

signal_power = max(cfg.target_amplitude)^2;

noise_power = signal_power / (10^(snr_db / 10));

noise = sqrt(noise_power / 2) * ...
    (randn(size(adc_cube)) + 1j * randn(size(adc_cube)));

adc_cube = adc_cube + noise;

end
