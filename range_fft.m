%% Hann window for Range FFT

n_range = (0:Ns-1)';

range_window = ...
    0.5 - 0.5*cos(2*pi*n_range/(Ns-1));


%% Range FFT without window

range_data_no_window = fft(adc_cube, [], 1);


%% Apply Range Hann window

adc_windowed = ...
    adc_cube .* reshape(range_window, Ns, 1, 1);

range_data = fft(adc_windowed, [], 1);


%% Positive range bins

half_N = floor(Ns / 2);

freq_axis = (0:half_N-1) * fs / Ns;

range_axis = c * freq_axis / (2 * S);


%% First chirp magnitude

range_mag_no_window = ...
    abs(range_data_no_window(1:half_N, 1, 1));

range_mag_window = ...
    abs(range_data(1:half_N, 1, 1));


%% Normalize for comparison

range_db_no_window = ...
    20*log10(range_mag_no_window / max(range_mag_no_window) + eps);

range_db_window = ...
    20*log10(range_mag_window / max(range_mag_window) + eps);


%% Plot before / after

figure;

plot(range_axis, range_db_no_window);

hold on;

plot(range_axis, range_db_window);

xlabel("Range [m]");
ylabel("Normalized Magnitude [dB]");

title("Range FFT - Hann Window Comparison");

legend("No Window", "Hann Window");

grid on;

xlim([20 65]);
ylim([-80 5]);