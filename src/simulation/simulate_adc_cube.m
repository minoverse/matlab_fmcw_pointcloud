%% Fast-time axis

n = 0:Ns-1;

t_fast = n / fs;


%% ADC cube

adc_cube = complex(zeros(Ns, Nrx, Nchirp));


%% Generate target signals

for m = 0:Nchirp-1

    for rx = 1:Nrx

        rx_signal = complex(zeros(1, Ns));

        for k = 1:Ntargets

            fast_phase = 2*pi * fb(k) * t_fast;

            slow_phase = 2*pi * fd(k) * m * Tchirp;


            azimuth = target_angle(k);

            elevation = target_elevation(k);


            direction_y = cosd(elevation) * sind(azimuth);

            direction_z = sind(elevation);


            spatial_phase = 2*pi/lambda * (rx_y(rx) * direction_y + rx_z(rx) * direction_z);


            target_signal = target_amplitude(k) * exp(1j * (fast_phase + slow_phase + spatial_phase));


            rx_signal = rx_signal + target_signal;

        end


        adc_cube(:, rx, m+1) = rx_signal.';

    end

end


%% Add complex Gaussian noise

snr_db = 20;

signal_power = max(target_amplitude)^2;

noise_power = signal_power / (10^(snr_db / 10));


noise = sqrt(noise_power / 2) * (randn(size(adc_cube)) + 1j * randn(size(adc_cube)));


adc_cube = adc_cube + noise;


%% Plot first chirp from RX1

adc_signal = adc_cube(:, 1, 1).';


figure;

plot(t_fast * 1e6, real(adc_signal));

hold on;

plot(t_fast * 1e6, imag(adc_signal));

xlabel("Time [us]");

ylabel("Amplitude");

legend("I", "Q");

title("First Chirp IQ Signal");

grid on;