%% Fast-time axis

n = 0:Ns-1;
t_fast = n / fs;

%% Target amplitude

A = 1;

%% ADC cube
% Dimension 1 = Fast Time samples
% Dimension 2 = RX channels
% Dimension 3 = Chirps / Slow Time

Nrx = 1;

adc_cube = complex(zeros(Ns, Nrx, Nchirp));

%% Generate beat signal for each chirp

for m = 0:Nchirp-1

    slow_phase = 2*pi * fd * m * Tchirp;

    adc_cube(:, 1, m+1) = ...
        A * exp(1j * (2*pi*fb*t_fast + slow_phase));

end

%% Plot first chirp I/Q

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