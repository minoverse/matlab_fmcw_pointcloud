
%% Group 3C - Noise-Only Monte Carlo Pfa

%% Experiment Settings

num_trials = 100;

threshold_values = [0, 0.5, 1, 3, 6];

num_thresholds = length(threshold_values);

noise_power = 0.01;

half_N = floor(Ns / 2);


%% CFAR Parameters

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

Pfa = 1e-5;

range_margin = Tr + Gr;
doppler_margin = Td + Gd;

range_width = 2 * range_margin + 1;
doppler_width = 2 * doppler_margin + 1;

total_cells = range_width * doppler_width;

guard_cells = (2 * Gr + 1) * (2 * Gd + 1);

Ntraining = total_cells - guard_cells;

alpha = Ntraining * (Pfa^(-1/Ntraining) - 1);


%% Valid CFAR CUTs

valid_range = range_margin+1:half_N-range_margin;

valid_doppler = doppler_margin+1:Nchirp-doppler_margin;

dense_cfar_cuts = length(valid_range) * length(valid_doppler);

fair_dense_ffts = half_N * Nrx;

assert(Ntraining == 216);
assert(dense_cfar_cuts == 25920);
assert(fair_dense_ffts == 2000);


%% Range and Doppler Windows

n_range = (0:Ns-1)';

range_window = 0.5 - 0.5*cos(2*pi*n_range/(Ns-1));

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);


%% CFAR Training Kernel

cfar_kernel = ones(range_width, doppler_width);

center_r = range_margin + 1;
center_d = doppler_margin + 1;

cfar_kernel(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = 0;

assert(sum(cfar_kernel(:)) == Ntraining);


%% Counters

dense_false_total = 0;

sparse_false_total = zeros(1, num_thresholds);

sparse_cut_total = zeros(1, num_thresholds);

roi_bin_total = zeros(1, num_thresholds);

expanded_bin_total = zeros(1, num_thresholds);

sparse_fft_total = zeros(1, num_thresholds);

empty_roi_trials = zeros(1, num_thresholds);

cfar_spotcheck_pass = false;

sparse_regression_pass = false;


%% Monte Carlo

for trial = 1:num_trials

    %% Generate Noise Only

    adc_cube = sqrt(noise_power / 2) * (randn(Ns,Nrx,Nchirp) + 1j*randn(Ns,Nrx,Nchirp));


    %% Range FFT

    adc_windowed = adc_cube .* reshape(range_window, Ns, 1, 1);

    range_data = fft(adc_windowed, [], 1);

    range_positive = range_data(1:half_N, :, :);


    %% Fair Dense Doppler FFT

    doppler_input = range_positive .* reshape(doppler_window, 1, 1, Nchirp);

    dense_doppler_cube = fft(doppler_input, Nchirp, 3);

    dense_doppler_cube = fftshift(dense_doppler_cube, 3);

    RD_power = reshape(sum(abs(dense_doppler_cube).^2, 2), half_N, Nchirp);


    %% Dense CFAR

    training_sum = conv2(RD_power, cfar_kernel, "same");

    dense_noise_map = zeros(half_N, Nchirp);

    dense_noise_map(valid_range, valid_doppler) = training_sum(valid_range, valid_doppler) / Ntraining;

    cfar_map = false(half_N, Nchirp);

    cfar_map(valid_range, valid_doppler) = RD_power(valid_range, valid_doppler) > alpha * dense_noise_map(valid_range, valid_doppler);

    dense_false_count = nnz(cfar_map);

    dense_false_total = dense_false_total + dense_false_count;


    %% Verify CFAR Calculation Against Original Window Formula

    if trial == 1

        check_cells = [11, 6; 100, 20; 250, 32; 490, 59];

        for j = 1:size(check_cells, 1)

            r = check_cells(j, 1);
            d = check_cells(j, 2);

            window = RD_power(r-range_margin:r+range_margin, d-doppler_margin:d+doppler_margin);

            training_cells = window(logical(cfar_kernel));

            reference_noise = mean(training_cells);

            reference_detection = RD_power(r,d) > alpha * reference_noise;

            noise_error = abs(reference_noise - dense_noise_map(r,d));

            noise_tolerance = 1e-10 * max(1, abs(reference_noise));

            assert(noise_error <= noise_tolerance);
            assert(reference_detection == cfar_map(r,d));

        end

        cfar_spotcheck_pass = true;

    end


    %% Threshold Sweep

    previous_roi_count = half_N + 1;

    previous_expanded_count = half_N + 1;

    previous_false_count = dense_false_count + 1;

    previous_cut_count = dense_cfar_cuts + 1;

    for t = 1:num_thresholds

        threshold_db = threshold_values(t);

        [roi_mask, roi_bins, expanded_bins] = roi_threshold_gate(range_positive, threshold_db, Tr, Gr);


        %% Valid ROI CFAR CUTs

        valid_roi_bins = roi_bins(roi_bins > range_margin & roi_bins <= half_N-range_margin);

        roi_cfar_cuts = length(valid_roi_bins) * length(valid_doppler);


        %% Sparse False Alarms Restricted to ROI

        sparse_false_count = nnz(cfar_map(valid_roi_bins, valid_doppler));


        %% Projected Sparse Doppler Workload

        sparse_fft_count = length(expanded_bins) * Nrx;


        %% Counter Validation

        assert(length(roi_bins) <= previous_roi_count);
        assert(length(expanded_bins) <= previous_expanded_count);

        assert(sparse_false_count <= previous_false_count);
        assert(roi_cfar_cuts <= previous_cut_count);

        assert(sparse_false_count <= dense_false_count);
        assert(sparse_false_count <= roi_cfar_cuts);

        assert(sparse_fft_count <= fair_dense_ffts);
        assert(roi_cfar_cuts <= dense_cfar_cuts);

        previous_roi_count = length(roi_bins);
        previous_expanded_count = length(expanded_bins);

        previous_false_count = sparse_false_count;
        previous_cut_count = roi_cfar_cuts;


        %% Accumulate Results

        sparse_false_total(t) = sparse_false_total(t) + sparse_false_count;

        sparse_cut_total(t) = sparse_cut_total(t) + roi_cfar_cuts;

        roi_bin_total(t) = roi_bin_total(t) + length(roi_bins);

        expanded_bin_total(t) = expanded_bin_total(t) + length(expanded_bins);

        sparse_fft_total(t) = sparse_fft_total(t) + sparse_fft_count;

        if isempty(roi_bins)

            empty_roi_trials(t) = empty_roi_trials(t) + 1;

        end


        %% One Actual Sparse CFAR Regression

        if trial == 1 && t == 1 && ~isempty(roi_bins)

            range_data_for_doppler = range_data;

            apply_clutter_removal = false;

            run("roi_cfar_margin_selection.m");

            run("selective_doppler_processing.m");

            run("roi_only_cfar.m");


            %% Compare Doppler Power

            dense_selected_power = RD_power(active_range_bins, :);

            power_difference = max(abs(sparse_RD_power(:) - dense_selected_power(:)));

            power_tolerance = 1e-10 * max(1, max(abs(dense_selected_power(:))));

            assert(power_difference <= power_tolerance);


            %% Compare CFAR Detections

            dense_selected_cfar = cfar_map(roi_bins, :);

            assert(isequal(roi_cfar_local, dense_selected_cfar));


            %% Compare CFAR Noise

            dense_selected_noise = dense_noise_map(roi_bins, :);

            noise_difference = max(abs(roi_noise_local(:) - dense_selected_noise(:)));

            noise_tolerance = 1e-10 * max(1, max(abs(dense_selected_noise(:))));

            assert(noise_difference <= noise_tolerance);

            sparse_regression_pass = true;

        end

    end

end


%% Calculate Empirical Pfa

dense_cut_total = num_trials * dense_cfar_cuts;

dense_empirical_pfa = dense_false_total / dense_cut_total;

expected_dense_false_alarms = dense_cut_total * Pfa;

sparse_full_grid_pfa = sparse_false_total / dense_cut_total;

sparse_evaluated_cut_pfa = NaN(1, num_thresholds);

for t = 1:num_thresholds

    if sparse_cut_total(t) > 0

        sparse_evaluated_cut_pfa(t) = sparse_false_total(t) / sparse_cut_total(t);

    end

end


%% Average Workload

mean_roi_bins = roi_bin_total / num_trials;

mean_expanded_bins = expanded_bin_total / num_trials;

mean_sparse_ffts = sparse_fft_total / num_trials;

mean_sparse_cfar_cuts = sparse_cut_total / num_trials;

fft_reduction_pct = 100 * (1 - mean_sparse_ffts / fair_dense_ffts);

cfar_reduction_pct = 100 * (1 - mean_sparse_cfar_cuts / dense_cfar_cuts);


%% Final Assertions

assert(cfar_spotcheck_pass, "Dense CFAR spot-check missing");

assert(sparse_regression_pass, "Sparse CFAR regression missing");

assert(dense_false_total >= 0);
assert(dense_false_total <= dense_cut_total);

assert(all(sparse_false_total <= dense_false_total));

assert(all(sparse_cut_total <= dense_cut_total));

assert(all(fft_reduction_pct >= 0));
assert(all(fft_reduction_pct <= 100));


%% Results Table

results_data = zeros(num_thresholds, 15);

for t = 1:num_thresholds

    results_data(t,:) = [
        threshold_values(t), num_trials, dense_false_total, dense_empirical_pfa, sparse_false_total(t), sparse_full_grid_pfa(t), sparse_evaluated_cut_pfa(t), mean_roi_bins(t), mean_expanded_bins(t), mean_sparse_ffts(t), mean_sparse_cfar_cuts(t), fft_reduction_pct(t), cfar_reduction_pct(t), empty_roi_trials(t), sparse_cut_total(t)
    ];

end

column_names = [
    "Threshold_dB", "Trials", "DenseFalseAlarms", "DenseEmpiricalPfa", "SparseFalseAlarms", "SparseFullGridPfa", "SparseEvaluatedCutPfa", "MeanROIBins", "MeanExpandedBins", "MeanProjectedSparseFFTs", "MeanSparseCFARCUTs", "ProjectedFFTReduction_pct", "CFARCutReduction_pct", "EmptyROITrials", "TotalSparseEvaluatedCUTs"
];

results_table = array2table(results_data, "VariableNames", column_names);


%% Output Directory

output_dir = "output/phase2/05_noise_only_monte_carlo_pfa";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

csv_file = fullfile(output_dir, "noise_only_pfa_threshold_results.csv");

writetable(results_table, csv_file);


%% Plot Results

fig = figure("Position", [100 100 1100 900]);

subplot(3,1,1);

plot(threshold_values, sparse_full_grid_pfa * 1e6, "-o", "LineWidth", 1.5);

hold on;

plot(threshold_values, dense_empirical_pfa * 1e6 * ones(1,num_thresholds), "--", "LineWidth", 1.5);

xlabel("ROI Threshold Above Median [dB]");
ylabel("False Alarms per Million CUTs");

title("Noise-Only False Alarm Rate - Full Grid Basis");

legend("Sparse", "Dense", "Location", "best");

grid on;


subplot(3,1,2);

plot(threshold_values, fft_reduction_pct, "-o", "LineWidth", 1.5);

xlabel("ROI Threshold Above Median [dB]");
ylabel("Doppler FFT Reduction [%]");

title("Projected Doppler FFT Workload Reduction Including CFAR Margin");

ylim([0 100]);
grid on;


subplot(3,1,3);

plot(threshold_values, mean_sparse_cfar_cuts / dense_cfar_cuts * 100, "-o", "LineWidth", 1.5);

xlabel("ROI Threshold Above Median [dB]");
ylabel("Evaluated CFAR CUTs [% of Dense]");

title("Sparse CFAR Evaluation Coverage");

ylim([0 100]);
grid on;


plot_file = fullfile(output_dir, "noise_only_pfa_vs_roi_threshold.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save TXT Summary

txt_file = fullfile(output_dir, "noise_only_monte_carlo_pfa_validation.txt");

file = fopen(txt_file, "w");

assert(file ~= -1, "Cannot create Group 3C TXT file");

fprintf(file, "=== Phase 2 Group 3C: Noise-Only Monte Carlo Pfa ===\n\n");

fprintf(file, "[Experiment]\n");

fprintf(file, "Trials = %d\n", num_trials);

fprintf(file, "Noise power = %.6f\n", noise_power);

fprintf(file, "Nominal CFAR Pfa = %.3e\n", Pfa);

fprintf(file, "Training cells = %d\n", Ntraining);

fprintf(file, "Dense valid CUTs per trial = %d\n", dense_cfar_cuts);

fprintf(file, "Total Dense CUT observations = %d\n", dense_cut_total);


fprintf(file, "\n[Dense CFAR]\n");

fprintf(file, "Observed false alarms = %d\n", dense_false_total);

fprintf(file, "Empirical Pfa = %.6e\n", dense_empirical_pfa);

fprintf(file, "Nominal-model expected false alarms = %.2f\n", expected_dense_false_alarms);


fprintf(file, "\n[Sparse ROI Threshold Sweep]\n");

for t = 1:num_thresholds

    fprintf(file, "\nThreshold = %.2f dB\n", threshold_values(t));

    fprintf(file, "Mean ROI bins = %.2f\n", mean_roi_bins(t));

    fprintf(file, "Mean expanded bins = %.2f\n", mean_expanded_bins(t));

    fprintf(file, "Mean projected Sparse FFTs = %.2f\n", mean_sparse_ffts(t));

    fprintf(file, "Mean CFAR CUTs = %.2f\n", mean_sparse_cfar_cuts(t));

    fprintf(file, "Sparse false alarms = %d\n", sparse_false_total(t));

    fprintf(file, "Sparse Pfa over full grid = %.6e\n", sparse_full_grid_pfa(t));

    fprintf(file, "Sparse Pfa over evaluated CUTs = %.6e\n", sparse_evaluated_cut_pfa(t));

    fprintf(file, "Empty ROI trials = %d\n", empty_roi_trials(t));

    fprintf(file, "Projected FFT workload reduction = %.2f%%\n", fft_reduction_pct(t));

    fprintf(file, "CFAR CUT reduction = %.2f%%\n", cfar_reduction_pct(t));

end


fprintf(file, "\n[Validation]\n");

fprintf(file, "Dense CFAR spot-check = PASS\n");

fprintf(file, "Sparse CFAR regression = PASS\n");

fprintf(file, "ROI monotonicity = PASS\n");

fprintf(file, "False alarm subset property = PASS\n");

fprintf(file, "Workload bounds = PASS\n");

fprintf(file, "Overall code validation = PASS\n");


fprintf(file, "\n[Limitations]\n");

fprintf(file, "All inputs are noise only. No target signals are present.\n");

fprintf(file, "Sparse full-grid Pfa includes skipped CUTs as non-detections.\n");

fprintf(file, "Sparse evaluated-CUT Pfa is NaN if no CUTs were evaluated.\n");

fprintf(file, "A zero sparse false-alarm count can result from an empty ROI.\n");

fprintf(file, "Sparse Monte Carlo detections use Dense CFAR restricted to ROI CUTs.\n");

fprintf(file, "An actual Sparse CFAR regression is checked on one trial.\n");

fprintf(file, "The CFAR scale assumes independent exponential power samples.\n");

fprintf(file, "Current 4-RX power summation and windowing can change empirical Pfa.\n");

fprintf(file, "Pfa is measured, not required to equal its nominal setting.\n");

fprintf(file, "Workload counters are not CPU latency or power measurements.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 3C Validation ---\n");

fprintf("Noise-only trials = %d\n", num_trials);

fprintf("Dense CFAR CUT observations = %d\n", dense_cut_total);

fprintf("Dense false alarms = %d\n", dense_false_total);

fprintf("Dense empirical Pfa = %.6e\n", dense_empirical_pfa);

fprintf("Nominal Pfa = %.6e\n", Pfa);

fprintf("Dense CFAR spot-check = PASS\n");

fprintf("Sparse CFAR regression = PASS\n");

fprintf("ROI monotonicity = PASS\n");

fprintf("False alarm subset property = PASS\n");

fprintf("Workload counter bounds = PASS\n");

fprintf("Group 3C code validation passed.\n");

fprintf("CSV saved: %s\n", csv_file);

fprintf("Plot saved: %s\n", plot_file);

fprintf("Summary saved: %s\n", txt_file);
