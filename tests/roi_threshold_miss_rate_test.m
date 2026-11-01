
%% Group 3B - ROI Threshold Sweep and Miss-Rate Test

%% Experiment Settings

threshold_values = [0, 3, 6, 9, 12];

target2_amplitudes = [0.7, 0.03, 0.005];

num_trials = 30;

range_tolerance = 1;
doppler_tolerance = 1;

num_scenes = length(target2_amplitudes);
num_thresholds = length(threshold_values);

assert(Ntargets == 2, "This test expects two targets");


%% CFAR Parameters

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

Pfa = 1e-5;

range_margin = Tr + Gr;
doppler_margin = Td + Gd;

total_cells = (2 * range_margin + 1) * (2 * doppler_margin + 1);
guard_cells = (2 * Gr + 1) * (2 * Gd + 1);

Ntraining = total_cells - guard_cells;

alpha = Ntraining * (Pfa^(-1 / Ntraining) - 1);

half_N = floor(Ns / 2);

valid_doppler_cuts = Nchirp - 2 * doppler_margin;

dense_cfar_cuts = (half_N - 2 * range_margin) * valid_doppler_cuts;

fair_dense_ffts = half_N * Nrx;

assert(Ntraining == 216);
assert(dense_cfar_cuts == 25920);
assert(fair_dense_ffts == 2000);


%% FFT Windows and Axes

n_range = (0:Ns-1)';

range_window = 0.5 - 0.5*cos(2*pi*n_range/(Ns-1));

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);

range_axis = c * (0:half_N-1) * fs / (2*S*Ns);

PRF = 1 / Tchirp;

fd_axis = (-Nchirp/2:Nchirp/2-1) * PRF / Nchirp;

velocity_axis = fd_axis * lambda / 2;


%% CFAR Training Mask

training_mask = true(2*range_margin+1, 2*doppler_margin+1);

center_r = range_margin + 1;
center_d = doppler_margin + 1;

training_mask(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = false;

assert(nnz(training_mask) == Ntraining);


%% Result Counters

roi_miss_count = zeros(num_scenes, num_thresholds, Ntargets);

dense_hit_count = zeros(num_scenes, Ntargets);

sparse_hit_count = zeros(num_scenes, num_thresholds, Ntargets);

dense_lost_count = zeros(num_scenes, num_thresholds, Ntargets);

roi_bin_sum = zeros(num_scenes, num_thresholds);

expanded_bin_sum = zeros(num_scenes, num_thresholds);

fft_count_sum = zeros(num_scenes, num_thresholds);

cfar_cut_sum = zeros(num_scenes, num_thresholds);

sparse_fft_regression_pass = false;


%% Monte Carlo Experiment

for scene = 1:num_scenes

    target_amplitude(2) = target2_amplitudes(scene);

    %% Generate Clean Target Signal Using Existing Simulator

    run("simulate_adc_cube.m");

    clean_adc_cube = adc_cube - noise;

    close(gcf);

    noise_power = max(target_amplitude)^2 / (10^(snr_db / 10));

    fprintf("\n--- Scene %d: Target 2 Amplitude = %.3f ---\n", scene, target_amplitude(2));


    for trial = 1:num_trials

        %% Independent Noise Realization

        trial_noise = sqrt(noise_power / 2) * (randn(size(clean_adc_cube)) + 1j*randn(size(clean_adc_cube)));

        adc_cube = clean_adc_cube + trial_noise;


        %% Dense Range FFT

        adc_windowed = adc_cube .* reshape(range_window, Ns, 1, 1);

        range_data = fft(adc_windowed, [], 1);

        range_positive = range_data(1:half_N, :, :);


        %% Fair Dense Doppler Reference

        range_doppler_windowed = range_positive .* reshape(doppler_window, 1, 1, Nchirp);

        dense_doppler_cube = fft(range_doppler_windowed, Nchirp, 3);

        dense_doppler_cube = fftshift(dense_doppler_cube, 3);

        RD_power = reshape(sum(abs(dense_doppler_cube).^2, 2), half_N, Nchirp);


        %% Local Dense CFAR Near Each Ground-Truth Target

        dense_target_ranges = cell(Ntargets, 1);
        target_range_bins = zeros(Ntargets, 1);

        for k = 1:Ntargets

            [~, r0] = min(abs(range_axis - target_range(k)));
            [~, d0] = min(abs(velocity_axis - target_velocity(k)));

            target_range_bins(k) = r0;

            r_first = max(r0-range_tolerance, range_margin+1);
            r_last = min(r0+range_tolerance, half_N-range_margin);

            d_first = max(d0-doppler_tolerance, doppler_margin+1);
            d_last = min(d0+doppler_tolerance, Nchirp-doppler_margin);

            detected_ranges = [];

            for r = r_first:r_last

                for d = d_first:d_last

                    window = RD_power(r-range_margin:r+range_margin, d-doppler_margin:d+doppler_margin);

                    training_cells = window(training_mask);

                    noise_estimate = mean(training_cells);

                    if RD_power(r, d) > alpha * noise_estimate

                        detected_ranges(end+1) = r;

                    end

                end

            end

            dense_target_ranges{k} = detected_ranges;

            if ~isempty(detected_ranges)

                dense_hit_count(scene, k) = dense_hit_count(scene, k) + 1;

            end

        end


        %% Threshold Sweep on Same Noise Realization

        previous_roi_count = half_N + 1;

        for t = 1:num_thresholds

            threshold_db = threshold_values(t);

            [roi_mask, roi_bins, expanded_bins] = roi_threshold_gate(range_positive, threshold_db, Tr, Gr);

            roi_count = length(roi_bins);

            expanded_count = length(expanded_bins);

            sparse_fft_count = expanded_count * Nrx;

            valid_roi_count = nnz(roi_mask(range_margin+1:half_N-range_margin));

            sparse_cfar_cuts = valid_roi_count * valid_doppler_cuts;


            %% Counter Validation

            assert(roi_count <= previous_roi_count);
            assert(expanded_count >= roi_count);
            assert(expanded_count <= half_N);
            assert(sparse_fft_count <= fair_dense_ffts);
            assert(sparse_cfar_cuts <= dense_cfar_cuts);

            previous_roi_count = roi_count;


            %% Accumulate Workload

            roi_bin_sum(scene, t) = roi_bin_sum(scene, t) + roi_count;

            expanded_bin_sum(scene, t) = expanded_bin_sum(scene, t) + expanded_count;

            fft_count_sum(scene, t) = fft_count_sum(scene, t) + sparse_fft_count;

            cfar_cut_sum(scene, t) = cfar_cut_sum(scene, t) + sparse_cfar_cuts;


            %% Target Coverage and Detection Retention

            for k = 1:Ntargets

                r0 = target_range_bins(k);

                nearby_bins = max(1, r0-range_tolerance):min(half_N, r0+range_tolerance);

                target_covered = any(roi_mask(nearby_bins));

                dense_hit = ~isempty(dense_target_ranges{k});

                sparse_hit = any(roi_mask(dense_target_ranges{k}));

                if ~target_covered

                    roi_miss_count(scene, t, k) = roi_miss_count(scene, t, k) + 1;

                end

                if sparse_hit

                    sparse_hit_count(scene, t, k) = sparse_hit_count(scene, t, k) + 1;

                end

                if dense_hit && ~sparse_hit

                    dense_lost_count(scene, t, k) = dense_lost_count(scene, t, k) + 1;

                end

                assert(~sparse_hit || dense_hit);

            end


            %% One Actual Selective FFT Regression

            if scene == 1 && trial == 1 && threshold_db == 6 && ~isempty(expanded_bins)

                range_data_for_doppler = range_positive;

                active_range_bins = expanded_bins;

                apply_clutter_removal = false;

                run("selective_doppler_processing.m");

                dense_selected = dense_doppler_cube(active_range_bins, :, :);

                fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_selected(:)));

                fft_tolerance = 1e-10 * max(1, max(abs(dense_selected(:))));

                assert(fft_difference <= fft_tolerance, "Selective FFT regression failed");

                sparse_fft_regression_pass = true;

            end

        end

    end

end


%% Calculate Rates

roi_miss_pct = 100 * roi_miss_count / num_trials;

dense_hit_pct = 100 * dense_hit_count / num_trials;

sparse_hit_pct = 100 * sparse_hit_count / num_trials;

dense_loss_pct = NaN(num_scenes, num_thresholds, Ntargets);

mean_roi_bins = roi_bin_sum / num_trials;

mean_expanded_bins = expanded_bin_sum / num_trials;

mean_fft_count = fft_count_sum / num_trials;

mean_cfar_cuts = cfar_cut_sum / num_trials;

fft_reduction_pct = 100 * (1 - mean_fft_count / fair_dense_ffts);

cfar_reduction_pct = 100 * (1 - mean_cfar_cuts / dense_cfar_cuts);


for scene = 1:num_scenes

    for t = 1:num_thresholds

        for k = 1:Ntargets

            if dense_hit_count(scene, k) > 0

                dense_loss_pct(scene, t, k) = 100 * dense_lost_count(scene, t, k) / dense_hit_count(scene, k);

            end

        end

    end

end


%% Statistical and Monotonicity Validation

assert(sparse_fft_regression_pass, "Selective FFT regression not executed");

assert(all(roi_miss_count(:) >= 0));
assert(all(roi_miss_count(:) <= num_trials));

assert(all(sparse_hit_count(:) >= 0));
assert(all(sparse_hit_count(:) <= num_trials));

for scene = 1:num_scenes

    for k = 1:Ntargets

        misses = squeeze(roi_miss_count(scene, :, k));

        hits = squeeze(sparse_hit_count(scene, :, k));

        assert(all(diff(misses) >= 0), "ROI misses must not decrease at higher thresholds");

        assert(all(diff(hits) <= 0), "Sparse target hits must not increase at higher thresholds");

        assert(all(hits <= dense_hit_count(scene, k)));

    end

end


%% Build Results Table

results_data = zeros(num_scenes * num_thresholds, 16);

row = 0;

for scene = 1:num_scenes

    for t = 1:num_thresholds

        row = row + 1;

        results_data(row, :) = [
            target2_amplitudes(scene), threshold_values(t), num_trials, roi_miss_pct(scene,t,1), roi_miss_pct(scene,t,2), dense_hit_pct(scene,1), dense_hit_pct(scene,2), sparse_hit_pct(scene,t,1), sparse_hit_pct(scene,t,2), dense_loss_pct(scene,t,2), mean_roi_bins(scene,t), mean_expanded_bins(scene,t), mean_fft_count(scene,t), mean_cfar_cuts(scene,t), fft_reduction_pct(scene,t), cfar_reduction_pct(scene,t)
        ];

    end

end

column_names = {
    "Target2Amplitude", "Threshold_dB", "Trials", "Target1ROIMiss_pct", "Target2ROIMiss_pct", "DenseTarget1Hit_pct", "DenseTarget2Hit_pct", "SparseTarget1Hit_pct", "SparseTarget2Hit_pct", "DenseConditionedT2Loss_pct", "MeanROIBins", "MeanExpandedBins", "MeanSparseFFTs", "MeanCFARCUTs", "DopplerReduction_pct", "CFARReduction_pct"
};

column_names = string(column_names);

results_table = array2table(results_data, "VariableNames", column_names);


%% Output Directory

output_dir = "output/phase2/04_roi_threshold_sweep_miss_rate";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

csv_file = fullfile(output_dir, "roi_threshold_sweep_miss_rate.csv");

writetable(results_table, csv_file);


%% Plot Results

fig = figure("Position", [100 100 1100 900]);

for scene = 1:num_scenes

    scene_name = sprintf("Target 2 amp = %.3f", target2_amplitudes(scene));


    subplot(3, 1, 1);

    hold on;

    plot(threshold_values, squeeze(roi_miss_pct(scene,:,2)), "-o", "DisplayName", scene_name);


    subplot(3, 1, 2);

    hold on;

    plot(threshold_values, squeeze(sparse_hit_pct(scene,:,2)), "-o", "DisplayName", "Sparse " + scene_name);

    plot(threshold_values, dense_hit_pct(scene,2) * ones(1,num_thresholds), "--", "DisplayName", "Dense " + scene_name);


    subplot(3, 1, 3);

    hold on;

    plot(threshold_values, fft_reduction_pct(scene,:), "-o", "DisplayName", scene_name);

end


subplot(3, 1, 1);

title("Target 2 - ROI Miss Rate");

xlabel("ROI Threshold Above Median [dB]");
ylabel("ROI Miss Rate [%]");

ylim([0 100]);
grid on;
legend("Location", "best");


subplot(3, 1, 2);

title("Target 2 - Local CFAR Target-Hit Rate");

xlabel("ROI Threshold Above Median [dB]");
ylabel("Target-Hit Rate [%]");

ylim([0 100]);
grid on;
legend("Location", "best");


subplot(3, 1, 3);

title("Projected Doppler FFT Workload Reduction (Including CFAR Margin)");

xlabel("ROI Threshold Above Median [dB]");
ylabel("Reduction vs Fair Dense [%]");

ylim([0 100]);
grid on;
legend("Location", "best");


plot_file = fullfile(output_dir, "roi_threshold_miss_rate_vs_doppler_workload.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save Text Summary

txt_file = fullfile(output_dir, "roi_threshold_sweep_validation.txt");

file = fopen(txt_file, "w");

assert(file ~= -1, "Cannot create Group 3B result file");

fprintf(file, "=== Phase 2 Group 3B: ROI Threshold and Miss-Rate Sweep ===\n\n");

fprintf(file, "Trials per amplitude = %d\n", num_trials);
fprintf(file, "Noise SNR = %.1f dB\n", snr_db);
fprintf(file, "Thresholds = 0, 3, 6, 9, 12 dB\n");
fprintf(file, "Target match tolerance = +/-1 range bin, +/-1 Doppler bin\n");
fprintf(file, "CFAR training cells = %d\n", Ntraining);
fprintf(file, "Fair Dense Doppler FFTs = %d\n", fair_dense_ffts);
fprintf(file, "Dense CFAR CUTs = %d\n", dense_cfar_cuts);

for scene = 1:num_scenes

    fprintf(file, "\n[Target 2 amplitude = %.3f]\n", target2_amplitudes(scene));

    fprintf(file, "Threshold | ROI Miss T2 | Dense Hit T2 | Sparse Hit T2 | Mean FFTs | FFT Reduction\n");

    for t = 1:num_thresholds

        fprintf(file, "%5.1f dB | %8.2f%% | %8.2f%% | %8.2f%% | %8.1f | %8.2f%%\n", threshold_values(t), roi_miss_pct(scene,t,2), dense_hit_pct(scene,2), sparse_hit_pct(scene,t,2), mean_fft_count(scene,t), fft_reduction_pct(scene,t));

    end

end

fprintf(file, "\n[Validation]\n");
fprintf(file, "Independent noise realizations = PASS\n");
fprintf(file, "Paired threshold comparison = PASS\n");
fprintf(file, "Threshold monotonicity = PASS\n");
fprintf(file, "Sparse target-hit subset of Dense = PASS\n");
fprintf(file, "Selective FFT spot-check = PASS\n");
fprintf(file, "Workload counter bounds = PASS\n");
fprintf(file, "Overall = PASS\n");

fprintf(file, "\n[Limitations]\n");
fprintf(file, "Target hit is evaluated on local CFAR CUTs near ground truth.\n");
fprintf(file, "Peak grouping and AoA are not part of the Monte Carlo hit metric.\n");
fprintf(file, "Sparse CFAR target hits use Dense CUT decisions restricted by ROI.\n");
fprintf(file, "Group 3A validated that ROI CFAR decisions equal Dense at the same CUTs.\n");
fprintf(file, "FFT and CFAR savings are workload projections, not latency or power.\n");
fprintf(file, "Thirty trials are exploratory and do not establish a low miss-rate guarantee.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 3B Validation ---\n");

fprintf("Trials per scene = %d\n", num_trials);

fprintf("Scene count = %d\n", num_scenes);

fprintf("Threshold monotonicity = PASS\n");

fprintf("Sparse target-hit subset of Dense = PASS\n");

fprintf("Selective FFT spot-check = PASS\n");

fprintf("Workload counter bounds = PASS\n");

fprintf("Group 3B validation passed.\n");

fprintf("CSV saved: %s\n", csv_file);
fprintf("Plot saved: %s\n", plot_file);
fprintf("Summary saved: %s\n", txt_file);
