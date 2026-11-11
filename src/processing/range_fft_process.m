
function [range_data, range_axis, range_data_no_window] = range_fft_process(adc_cube, cfg)

%% Range Hann window

n_range = (0:cfg.Ns-1)';

range_window = 0.5 - 0.5*cos(2*pi*n_range/(cfg.Ns-1));


%% Range FFT without window

range_data_no_window = fft(adc_cube, [], 1);


%% Range FFT with Hann window

adc_windowed = adc_cube .* reshape(range_window, cfg.Ns, 1, 1);

range_data = fft(adc_windowed, [], 1);


%% Positive range bins

half_N = floor(cfg.Ns / 2);

freq_axis = (0:half_N-1) * cfg.fs / cfg.Ns;

range_axis = cfg.c * freq_axis / (2 * cfg.S);

end
