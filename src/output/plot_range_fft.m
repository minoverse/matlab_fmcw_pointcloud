
function plot_range_fft(range_data, range_data_no_window, range_axis)

%% First chirp magnitude

half_N = length(range_axis);

range_mag_no_window = abs(range_data_no_window(1:half_N, 1, 1));

range_mag_window = abs(range_data(1:half_N, 1, 1));


%% Normalize for comparison

range_db_no_window = 20*log10(range_mag_no_window / max(range_mag_no_window) + eps);

range_db_window = 20*log10(range_mag_window / max(range_mag_window) + eps);


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

end
