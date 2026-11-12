function cfar_pfa_calibration_test()

cfg = radar_config();
c = cfg.c;
fc = cfg.fc;
B = cfg.B;
Tchirp = cfg.Tchirp;
fs = cfg.fs;
Nchirp = cfg.Nchirp;
target_range = cfg.target_range;
target_velocity = cfg.target_velocity;
target_angle = cfg.target_angle;
target_elevation = cfg.target_elevation;
target_amplitude = cfg.target_amplitude;
Ntargets = cfg.Ntargets;
S = cfg.S;
Ns = cfg.Ns;
range_resolution = cfg.range_resolution;
tau = cfg.tau;
fb = cfg.fb;
lambda = cfg.lambda;
fd = cfg.fd;
Nrx_az = cfg.Nrx_az;
Nrx_el = cfg.Nrx_el;
Nrx = cfg.Nrx;
d_az = cfg.d_az;
d_el = cfg.d_el;
rx_y = cfg.rx_y;
rx_z = cfg.rx_z;
NangleFFT_az = cfg.NangleFFT_az;
NangleFFT_el = cfg.NangleFFT_el;


%% Group 3C-2 - CFAR Pfa Calibration

%% Experiment Settings

num_ideal_trials = 200;
num_calibration_trials = 200;
num_validation_trials = 300;

noise_power = 0.01;

target_pfa = 1e-5;

half_N = floor(Ns / 2);

assert(Ns == 1000);
assert(Nrx == 4);
assert(Nchirp == 64);


%% CFAR Parameters

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

range_margin = Tr + Gr;
doppler_margin = Td + Gd;

Ntraining = (2*range_margin+1) * (2*doppler_margin+1) - (2*Gr+1) * (2*Gd+1);

alpha_original = Ntraining * (target_pfa^(-1/Ntraining) - 1);

valid_range = range_margin+1:half_N-range_margin;

valid_doppler = doppler_margin+1:Nchirp-doppler_margin;

dense_cfar_cuts = length(valid_range) * length(valid_doppler);

assert(Ntraining == 216);
assert(dense_cfar_cuts == 25920);


%% CFAR Training Kernel

cfar_kernel = ones(2*range_margin+1, 2*doppler_margin+1);

center_r = range_margin + 1;
center_d = doppler_margin + 1;

cfar_kernel(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = 0;

assert(sum(cfar_kernel(:)) == Ntraining);


%% FFT Windows

n_range = (0:Ns-1)';

range_window = 0.5 - 0.5*cos(2*pi*n_range/(Ns-1));

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);


%% Alpha Candidates

alpha_candidates = [3:0.1:8, alpha_original];

num_candidates = length(alpha_candidates);

calibration_false = zeros(1, num_candidates);


%% Test 1 - Ideal Independent Exponential Noise

ideal_false_total = 0;

for trial = 1:num_ideal_trials

    ideal_power = -log(rand(half_N, Nchirp));

    training_sum = conv2(ideal_power, cfar_kernel, "same");

    noise_estimate = training_sum(valid_range, valid_doppler) / Ntraining;

    cut_power = ideal_power(valid_range, valid_doppler);

    ideal_false_total = ideal_false_total + nnz(cut_power > alpha_original * noise_estimate);

end

ideal_total_cuts = num_ideal_trials * dense_cfar_cuts;

ideal_empirical_pfa = ideal_false_total / ideal_total_cuts;

ideal_expected_false = ideal_total_cuts * target_pfa;

ideal_check = ideal_empirical_pfa >= 0.5*target_pfa && ideal_empirical_pfa <= 1.5*target_pfa;


%% Test 2 - Actual 4-RX Calibration

real_1rx_false_total = 0;

for trial = 1:num_calibration_trials

    %% Generate Noise-Only ADC Cube

    adc_cube = sqrt(noise_power/2) * (randn(Ns,Nrx,Nchirp) + 1j*randn(Ns,Nrx,Nchirp));


    %% Range FFT

    adc_windowed = adc_cube .* reshape(range_window, Ns, 1, 1);

    range_data = fft(adc_windowed, [], 1);

    range_positive = range_data(1:half_N, :, :);


    %% Doppler FFT

    doppler_input = range_positive .* reshape(doppler_window, 1, 1, Nchirp);

    doppler_cube = fft(doppler_input, Nchirp, 3);

    doppler_cube = fftshift(doppler_cube, 3);


    %% 1-RX Processed Reference

    power_1rx = reshape(abs(doppler_cube(:,1,:)).^2, half_N, Nchirp);

    noise_sum_1rx = conv2(power_1rx, cfar_kernel, "same");

    noise_1rx = noise_sum_1rx(valid_range,valid_doppler) / Ntraining;

    ratio_1rx = power_1rx(valid_range,valid_doppler) ./ noise_1rx;

    real_1rx_false_total = real_1rx_false_total + nnz(ratio_1rx > alpha_original);


    %% Actual 4-RX Power

    power_4rx = reshape(sum(abs(doppler_cube).^2, 2), half_N, Nchirp);

    noise_sum_4rx = conv2(power_4rx, cfar_kernel, "same");

    noise_4rx = noise_sum_4rx(valid_range,valid_doppler) / Ntraining;

    ratio_4rx = power_4rx(valid_range,valid_doppler) ./ noise_4rx;


    %% CFAR Formula Spot-Check

    if trial == 1

        check_r = [11, 100, 250, 490];
        check_d = [6, 20, 32, 59];

        for j = 1:length(check_r)

            r = check_r(j);
            d = check_d(j);

            window = power_4rx(r-range_margin:r+range_margin, d-doppler_margin:d+doppler_margin);

            training_cells = window(logical(cfar_kernel));

            direct_noise = mean(training_cells);

            convolution_noise = noise_sum_4rx(r,d) / Ntraining;

            tolerance = 1e-10 * max(1, abs(direct_noise));

            assert(abs(direct_noise - convolution_noise) <= tolerance);

        end

    end


    %% Alpha Sweep

    for a = 1:num_candidates

        calibration_false(a) = calibration_false(a) + nnz(ratio_4rx > alpha_candidates(a));

    end

end


%% Calculate Calibration Pfa

calibration_total_cuts = num_calibration_trials * dense_cfar_cuts;

calibration_pfa = calibration_false / calibration_total_cuts;

real_1rx_empirical_pfa = real_1rx_false_total / calibration_total_cuts;


%% Select Best Alpha Using Calibration Data Only

calibration_error = abs(calibration_pfa - target_pfa);

[~, best_index] = min(calibration_error);

alpha_calibrated = alpha_candidates(best_index);

selected_calibration_pfa = calibration_pfa(best_index);

selected_calibration_false = calibration_false(best_index);


%% Test 3 - Independent 4-RX Validation

validation_original_false = 0;

validation_calibrated_false = 0;

for trial = 1:num_validation_trials

    %% New Independent Noise

    adc_cube = sqrt(noise_power/2) * (randn(Ns,Nrx,Nchirp) + 1j*randn(Ns,Nrx,Nchirp));


    %% Range FFT

    adc_windowed = adc_cube .* reshape(range_window, Ns, 1, 1);

    range_data = fft(adc_windowed, [], 1);

    range_positive = range_data(1:half_N, :, :);


    %% Doppler FFT

    doppler_input = range_positive .* reshape(doppler_window, 1, 1, Nchirp);

    doppler_cube = fft(doppler_input, Nchirp, 3);

    doppler_cube = fftshift(doppler_cube, 3);


    %% 4-RX CFAR Ratio

    power_4rx = reshape(sum(abs(doppler_cube).^2, 2), half_N, Nchirp);

    noise_sum_4rx = conv2(power_4rx, cfar_kernel, "same");

    noise_4rx = noise_sum_4rx(valid_range,valid_doppler) / Ntraining;

    ratio_4rx = power_4rx(valid_range,valid_doppler) ./ noise_4rx;


    %% Compare Original and Calibrated Thresholds

    validation_original_false = validation_original_false + nnz(ratio_4rx > alpha_original);

    validation_calibrated_false = validation_calibrated_false + nnz(ratio_4rx > alpha_calibrated);

end


%% Independent Validation Results

validation_total_cuts = num_validation_trials * dense_cfar_cuts;

validation_original_pfa = validation_original_false / validation_total_cuts;

validation_calibrated_pfa = validation_calibrated_false / validation_total_cuts;

validation_expected_false = validation_total_cuts * target_pfa;


%% Calibration Acceptance

pfa_lower_limit = 0.5 * target_pfa;
pfa_upper_limit = 1.5 * target_pfa;

calibration_pass = validation_calibrated_pfa >= pfa_lower_limit && validation_calibrated_pfa <= pfa_upper_limit;

if ideal_check
    ideal_status = "PASS";
else
    ideal_status = "NEEDS REVIEW";
end

if calibration_pass
    calibration_status = "PASS";
else
    calibration_status = "NEEDS REVIEW";
end


%% Counter Validation

assert(all(diff(calibration_false) <= 0));

assert(alpha_calibrated > 0);

assert(validation_calibrated_false >= 0);

assert(validation_calibrated_false <= validation_total_cuts);

assert(validation_original_false <= validation_total_cuts);


%% Save Directory

output_dir = "output/phase2/06_cfar_pfa_calibration";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Save CSV

results_data = [alpha_candidates(:), calibration_false(:), calibration_pfa(:), calibration_error(:)];

column_names = {'Alpha', 'CalibrationFalseAlarms', 'CalibrationEmpiricalPfa', 'AbsolutePfaError'};

results_table = array2table(results_data, "VariableNames", column_names);

csv_file = fullfile(output_dir, "cfar_alpha_sweep_calibration.csv");

writetable(results_table, csv_file);


%% Plot

fig = figure("Position", [100 100 1100 750]);

subplot(2,1,1);

plot_floor = 0.5 / calibration_total_cuts;

semilogy(alpha_candidates, max(calibration_pfa, plot_floor), "-o");

hold on;

yline(target_pfa, "r--");

xline(alpha_calibrated, "g--");

xlabel("CFAR Alpha");
ylabel("Empirical Pfa");

title("4-RX CFAR Alpha Calibration");

legend("Calibration Pfa", "Nominal Pfa", "Selected Alpha");

grid on;


subplot(2,1,2);

pfa_comparison = [target_pfa, validation_original_pfa, validation_calibrated_pfa] * 1e6;

bar(pfa_comparison);

set(gca, "XTickLabel", {"Nominal", "Original 4-RX", "Calibrated 4-RX"});

ylabel("False Alarms per Million CUTs");

title("Independent Validation - 300 Trials");

grid on;


plot_file = fullfile(output_dir, "cfar_pfa_calibration_and_holdout_validation.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save TXT Results

txt_file = fullfile(output_dir, "cfar_pfa_calibration_validation.txt");

file = fopen(txt_file, "w");

assert(file ~= -1, "Cannot create calibration result file");

fprintf(file, "=== Phase 2 Group 3C-2: CFAR Pfa Calibration ===\n\n");

fprintf(file, "[CFAR Settings]\n");
fprintf(file, "Nominal Pfa = %.6e\n", target_pfa);
fprintf(file, "Training cells = %d\n", Ntraining);
fprintf(file, "Valid CUTs per trial = %d\n", dense_cfar_cuts);
fprintf(file, "Original alpha = %.6f\n", alpha_original);

fprintf(file, "\n[Ideal Independent 1-RX Test]\n");
fprintf(file, "Trials = %d\n", num_ideal_trials);
fprintf(file, "Expected false alarms = %.2f\n", ideal_expected_false);
fprintf(file, "Observed false alarms = %d\n", ideal_false_total);
fprintf(file, "Empirical Pfa = %.6e\n", ideal_empirical_pfa);
fprintf(file, "Ideal formula check = %s\n", ideal_status);

fprintf(file, "\n[Actual FFT 1-RX Diagnostic]\n");
fprintf(file, "Trials = %d\n", num_calibration_trials);
fprintf(file, "Observed false alarms = %d\n", real_1rx_false_total);
fprintf(file, "Empirical Pfa = %.6e\n", real_1rx_empirical_pfa);

fprintf(file, "\n[4-RX Calibration]\n");
fprintf(file, "Calibration trials = %d\n", num_calibration_trials);
fprintf(file, "Selected alpha = %.6f\n", alpha_calibrated);
fprintf(file, "Calibration false alarms = %d\n", selected_calibration_false);
fprintf(file, "Calibration empirical Pfa = %.6e\n", selected_calibration_pfa);

fprintf(file, "\n[Independent 4-RX Validation]\n");
fprintf(file, "Validation trials = %d\n", num_validation_trials);
fprintf(file, "Total CUT observations = %d\n", validation_total_cuts);
fprintf(file, "Nominal expected false alarms = %.2f\n", validation_expected_false);
fprintf(file, "Original alpha false alarms = %d\n", validation_original_false);
fprintf(file, "Original alpha empirical Pfa = %.6e\n", validation_original_pfa);
fprintf(file, "Calibrated alpha false alarms = %d\n", validation_calibrated_false);
fprintf(file, "Calibrated empirical Pfa = %.6e\n", validation_calibrated_pfa);
fprintf(file, "Acceptance range = %.6e to %.6e\n", pfa_lower_limit, pfa_upper_limit);
fprintf(file, "Calibration validation = %s\n", calibration_status);

fprintf(file, "\n[Notes]\n");
fprintf(file, "Ideal test uses independent exponential power samples.\n");
fprintf(file, "Actual 4-RX test uses Range FFT, Doppler FFT, Hann windows and RX power summation.\n");
fprintf(file, "Calibration and validation use separate noise realizations.\n");
fprintf(file, "Acceptance range is an exploratory engineering tolerance, not a formal confidence interval.\n");
fprintf(file, "The original CFAR implementation is not modified by this test.\n");

fclose(file);


%% Final Output

fprintf("\n--- Group 3C-2 CFAR Calibration ---\n");

fprintf("Original alpha = %.6f\n", alpha_original);

fprintf("Ideal 1-RX empirical Pfa = %.6e\n", ideal_empirical_pfa);

fprintf("Ideal formula check = %s\n", ideal_status);

fprintf("Actual FFT 1-RX empirical Pfa = %.6e\n", real_1rx_empirical_pfa);

fprintf("Calibrated 4-RX alpha = %.6f\n", alpha_calibrated);

fprintf("Calibration empirical Pfa = %.6e\n", selected_calibration_pfa);

fprintf("Independent validation empirical Pfa = %.6e\n", validation_calibrated_pfa);

fprintf("Original 4-RX validation empirical Pfa = %.6e\n", validation_original_pfa);

fprintf("Calibration validation = %s\n", calibration_status);

fprintf("CSV saved: %s\n", csv_file);
fprintf("Plot saved: %s\n", plot_file);
fprintf("Summary saved: %s\n", txt_file);

end
