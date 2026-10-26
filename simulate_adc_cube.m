%% Fast-time axis

n = 0:Ns-1;
t_fast = n / fs;


%% ADC cube
% Dimension 1 = Fast Time samples
% Dimension 2 = RX channels
% Dimension 3 = Chirps / Slow Time

Nrx = 1;

adc_cube = complex(zeros(Ns, Nrx, Nchirp));


%% Generate multiple target beat signals

for m = 0:Nchirp-1

    chirp_signal = complex(zeros(1, Ns));

    for k = 1:Ntargets

        fast_phase = 2*pi * fb(k) * t_fast;

        slow_phase = 2*pi * fd(k) * m * Tchirp;

        target_signal = target_amplitude(k) * ...
            exp(1j * (fast_phase + slow_phase));

        chirp_signal = chirp_signal + target_signal;

    end

    adc_cube(:, 1, m+1) = chirp_signal.';

end


%% Add complex Gaussian noise

snr_db = 20;

signal_power = max(target_amplitude)^2;

noise_power = signal_power / (10^(snr_db / 10));

noise = sqrt(noise_power / 2) * ...
    (randn(size(adc_cube)) + 1j * randn(size(adc_cube)));

adc_cube = adc_cube + noise;


%% Plot first chirp I/Q

adc_signal = adc_cube(:, 1, 1).';

figure;

plot(t_fast * 1e6, real(adc_signal));

hold on;

plot(t_fast * 1e6, imag(adc_signal));

xlabel("Time [us]");
ylabel("Amplitude");

legend("I", "Q");

title("First Chirp IQ Signal - Multiple Targets");

grid on;