
function plot_first_chirp_iq(adc_cube, t_fast)

%% First chirp from RX1

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

end
