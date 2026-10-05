% 1. FFT
range_spectrum = fft(adc_signal);

% 2. Magnitude
range_mag = abs(range_spectrum);

% 3. Frequency bins
freq_axis = (0:Ns-1) * fs / Ns;

% 4. Frequency → Range
range_axis = c * freq_axis / (2 * S);

%% Range FFT

range_spectrum = fft(adc_signal);

% Use only positive frequencies
half_N = floor(Ns / 2);

range_spectrum_pos = range_spectrum(1:half_N);
range_mag = abs(range_spectrum_pos);

%% Frequency axis

freq_axis = (0:half_N-1) * fs / Ns;

%% Frequency -> Range

range_axis = c * freq_axis / (2 * S);

%% Find peak

[peak_mag, peak_idx] = max(range_mag);

estimated_fb = freq_axis(peak_idx);
estimated_range = range_axis(peak_idx);

range_error = estimated_range - target_range;

%% Print result

fprintf("\n--- Range FFT Result ---\n");
fprintf("Ground Truth Range = %.3f m\n", target_range);
fprintf("Estimated Beat     = %.3f MHz\n", estimated_fb / 1e6);
fprintf("Estimated Range    = %.3f m\n", estimated_range);
fprintf("Range Error        = %.3f m\n", range_error);

%% Plot

figure;

plot(range_axis, range_mag);

xlabel("Range [m]");
ylabel("Magnitude");
title("Range FFT");

grid on;

xlim([0 60]);