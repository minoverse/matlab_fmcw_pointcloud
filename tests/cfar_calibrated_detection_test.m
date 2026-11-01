
%% Group 3D - Calibrated CFAR Detection Validation

%% Experiment Settings

num_trials = 200;

target2_amplitudes = [0.7, 0.03, 0.005];

threshold_values = [3, 6];

alpha_calibrated = 4.8;

target_pfa = 1e-5;

range_tolerance = 1;
doppler_tolerance = 1;

num_scenes = length(target2_amplitudes);
num_thresholds = length(threshold_values);
num_alpha = 2;

assert(Ntargets == 2);


%% CFAR Parameters

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

range_margin = Tr + Gr;
doppler_margin = Td + Gd;

Ntraining = (2*range_margin+1)*(2*doppler_margin+1) - (2*Gr+1)*(2*Gd+1);

alpha_original = Ntraining * (target_pfa^(-1/Ntraining) - 1);

alpha_values = [alpha_original, alpha_calibrated];

half_N = floor(Ns/2);

valid_range = range_margin+1:half_N-range_margin;
valid_doppler = doppler_margin+1:Nchirp-doppler_margin;

dense_cfar_cuts = length(valid_range) * length(valid_doppler);
fair_dense_ffts = half_N * Nrx;

assert(Ntraining == 216);
assert(dense_cfar_cuts == 25920);
assert(fair_dense_ffts == 2000);
assert(alpha_calibrated < alpha_original);


%% CFAR Training Kernel

cfar_kernel = ones(2*range_margin+1, 2*doppler_margin+1);

center_r = range_margin + 1;
center_d = doppler_margin + 1;

cfar_kernel(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = 0;

training_mask = logical(cfar_kernel);

assert(nnz(training_mask) == Ntraining);


%% FFT Windows

n_range = (0:Ns-1)';

range_window = 0.5 - 0.5*cos(2*pi*n_range/(Ns-1));

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5 - 0.5*cos(2*pi*n_doppler/Nchirp);


%% Axes

range_axis = c * (0:half_N-1) * fs / (2*S*Ns);

PRF = 1/Tchirp;

velocity_axis = (-Nchirp/2:Nchirp/2-1) * PRF * lambda / (2*Nchirp);


%% Result Counters

dense_hit_count = zeros(num_scenes, num_alpha, Ntargets);

sparse_hit_count = zeros(num_scenes, num_thresholds, num_alpha, Ntargets);

roi_coverage_count = zeros(num_scenes, num_thresholds, Ntargets);

dense_offtarget_count = zeros(num_scenes, num_alpha);

sparse_offtarget_count = zeros(num_scenes, num_thresholds, num_alpha);

roi_bin_sum = zeros(num_scenes, num_thresholds);

expanded_bin_sum = zeros(num_scenes, num_thresholds);

fft_count_sum = zeros(num_scenes, num_thresholds);

cfar_cut_sum = zeros(num_scenes, num_thresholds);

sparse_regression_pass = false;


%% Monte Carlo

for scene = 1:num_scenes

    target_amplitude(2) = target2_amplitudes(scene);


    %% Generate Clean Target Signal

    run("simulate_adc_cube.m");

    clean_adc_cube = adc_cube - noise;

    close(gcf);


    %% Ground Truth Matching Regions

    target_range_bins = zeros(Ntargets, 1);

    target_doppler_bins = zeros(Ntargets, 1);

    target_match_mask = false(half_N, Nchirp);

    for k = 1:Ntargets

        [~, r0] = min(abs(range_axis - target_range(k)));

        [~, d0] = min(abs(velocity_axis - target_velocity(k)));

        target_range_bins(k) = r0;
        target_doppler_bins(k) = d0;

        r_first = max(range_margin+1, r0-range_tolerance);
        r_last = min(half_N-range_margin, r0+range_tolerance);

        d_first = max(doppler_margin+1, d0-doppler_tolerance);
        d_last = min(Nchirp-doppler_margin, d0+doppler_tolerance);

        target_match_mask(r_first:r_last, d_first:d_last) = true;

    end

    noise_power = max(target_amplitude)^2 / 10^(snr_db/10);

    fprintf("\n--- Scene %d: Target 2 amplitude = %.3f ---\n", scene, target_amplitude(2));


    for trial = 1:num_trials

        %% Independent Noise

        trial_noise = sqrt(noise_power/2) * (randn(size(clean_adc_cube)) + 1j*randn(size(clean_adc_cube)));

        adc_cube = clean_adc_cube + trial_noise;


        %% Range FFT

        adc_windowed = adc_cube .* reshape(range_window, Ns, 1, 1);

        range_data = fft(adc_windowed, [], 1);

        range_positive = range_data(1:half_N, :, :);


        %% Fair Dense Doppler FFT

        doppler_input = range_positive .* reshape(doppler_window, 1, 1, Nchirp);

        dense_doppler_cube = fft(doppler_input, Nchirp, 3);

        dense_doppler_cube = fftshift(dense_doppler_cube, 3);

        RD_power = reshape(sum(abs(dense_doppler_cube).^2, 2), half_N, Nchirp);


        %% CFAR Noise Estimation

        training_sum = conv2(RD_power, cfar_kernel, "same");

        noise_map = training_sum / Ntraining;


        %% Original and Calibrated Dense CFAR

        dense_maps = false(half_N, Nchirp, num_alpha);

        for a = 1:num_alpha

            dense_maps(valid_range, valid_doppler, a) = RD_power(valid_range, valid_doppler) > alpha_values(a) * noise_map(valid_range, valid_doppler);

            dense_map = dense_maps(:,:,a);

            dense_offtarget_count(scene,a) = dense_offtarget_count(scene,a) + nnz(dense_map & ~target_match_mask);

            for k = 1:Ntargets

                r0 = target_range_bins(k);
                d0 = target_doppler_bins(k);

                r_near = max(1,r0-range_tolerance):min(half_N,r0+range_tolerance);

                d_near = max(1,d0-doppler_tolerance):min(Nchirp,d0+doppler_tolerance);

                if any(any(dense_map(r_near,d_near)))

                    dense_hit_count(scene,a,k) = dense_hit_count(scene,a,k) + 1;

                end

            end

        end


        %% ROI Threshold Sweep

        previous_roi_count = half_N + 1;

        for t = 1:num_thresholds

            [roi_mask, roi_bins, expanded_bins] = roi_threshold_gate(range_positive, threshold_values(t), Tr, Gr);

            assert(length(roi_bins) <= previous_roi_count);

            previous_roi_count = length(roi_bins);


            %% Workload Counters

            roi_count = length(roi_bins);

            expanded_count = length(expanded_bins);

            sparse_fft_count = expanded_count * Nrx;

            valid_roi_count = nnz(roi_mask(valid_range));

            sparse_cfar_cuts = valid_roi_count * length(valid_doppler);

            assert(sparse_fft_count <= fair_dense_ffts);
            assert(sparse_cfar_cuts <= dense_cfar_cuts);

            roi_bin_sum(scene,t) = roi_bin_sum(scene,t) + roi_count;

            expanded_bin_sum(scene,t) = expanded_bin_sum(scene,t) + expanded_count;

            fft_count_sum(scene,t) = fft_count_sum(scene,t) + sparse_fft_count;

            cfar_cut_sum(scene,t) = cfar_cut_sum(scene,t) + sparse_cfar_cuts;


            %% Target ROI Coverage

            for k = 1:Ntargets

                r0 = target_range_bins(k);

                r_near = max(1,r0-range_tolerance):min(half_N,r0+range_tolerance);

                if any(roi_mask(r_near))

                    roi_coverage_count(scene,t,k) = roi_coverage_count(scene,t,k) + 1;

                end

            end


            %% Sparse CFAR Based on Same ROI CUTs

            sparse_mask = repmat(roi_mask, 1, Nchirp);

            for a = 1:num_alpha

                dense_map = dense_maps(:,:,a);

                sparse_map = dense_map & sparse_mask;

                assert(nnz(sparse_map & ~dense_map) == 0);

                sparse_offtarget_count(scene,t,a) = sparse_offtarget_count(scene,t,a) + nnz(sparse_map & ~target_match_mask);

                for k = 1:Ntargets

                    r0 = target_range_bins(k);
                    d0 = target_doppler_bins(k);

                    r_near = max(1,r0-range_tolerance):min(half_N,r0+range_tolerance);

                    d_near = max(1,d0-doppler_tolerance):min(Nchirp,d0+doppler_tolerance);

                    if any(any(sparse_map(r_near,d_near)))

                        sparse_hit_count(scene,t,a,k) = sparse_hit_count(scene,t,a,k) + 1;

                    end

                end

            end


            %% Actual Sparse FFT and CFAR Spot-Check

            if scene == 1 && trial == 1 && t == 1

                range_data_for_doppler = range_data;

                apply_clutter_removal = false;

                active_range_bins = expanded_bins;

                run("selective_doppler_processing.m");

                dense_selected = dense_doppler_cube(active_range_bins,:,:);

                fft_difference = max(abs(sparse_range_doppler_rx(:) - dense_selected(:)));

                assert(fft_difference <= 1e-10 * max(1,max(abs(dense_selected(:)))));


                %% Compare Sparse CFAR for Both Alpha Values

                range_lookup = zeros(half_N,1);

                range_lookup(active_range_bins) = 1:length(active_range_bins);

                for i = 1:length(roi_bins)

                    r = roi_bins(i);

                    if r <= range_margin || r > half_N-range_margin
                        continue;
                    end

                    local_rows = range_lookup(r-range_margin:r+range_margin);

                    assert(all(local_rows > 0));

                    for d = valid_doppler

                        window = sparse_RD_power(local_rows, d-doppler_margin:d+doppler_margin);

                        sparse_noise = mean(window(training_mask));

                        sparse_cut = sparse_RD_power(range_lookup(r),d);

                        for a = 1:num_alpha

                            sparse_detection = sparse_cut > alpha_values(a)*sparse_noise;

                            assert(sparse_detection == dense_maps(r,d,a));

                        end

                    end

                end

                sparse_regression_pass = true;

            end

        end

    end

end


%% Validation

assert(sparse_regression_pass);

for scene = 1:num_scenes

    for k = 1:Ntargets

        assert(dense_hit_count(scene,2,k) >= dense_hit_count(scene,1,k));

        for t = 1:num_thresholds

            assert(sparse_hit_count(scene,t,2,k) >= sparse_hit_count(scene,t,1,k));

            for a = 1:num_alpha

                assert(sparse_hit_count(scene,t,a,k) <= dense_hit_count(scene,a,k));

            end

        end

    end

end


%% Prepare Results

results_data = zeros(num_scenes*num_thresholds*num_alpha, 17);

row = 0;

for scene = 1:num_scenes

    for t = 1:num_thresholds

        for a = 1:num_alpha

            row = row + 1;

            if dense_hit_count(scene,a,2) > 0

                conditional_loss = 100 * (dense_hit_count(scene,a,2) - sparse_hit_count(scene,t,a,2)) / dense_hit_count(scene,a,2);

            else

                conditional_loss = NaN;

            end

            results_data(row,:) = [target2_amplitudes(scene), threshold_values(t), alpha_values(a), num_trials, 100*dense_hit_count(scene,a,1)/num_trials, 100*dense_hit_count(scene,a,2)/num_trials, 100*sparse_hit_count(scene,t,a,1)/num_trials, 100*sparse_hit_count(scene,t,a,2)/num_trials, 100*roi_coverage_count(scene,t,2)/num_trials, conditional_loss, roi_bin_sum(scene,t)/num_trials, expanded_bin_sum(scene,t)/num_trials, fft_count_sum(scene,t)/num_trials, cfar_cut_sum(scene,t)/num_trials, 100*(1-fft_count_sum(scene,t)/(num_trials*fair_dense_ffts)), dense_offtarget_count(scene,a)/num_trials, sparse_offtarget_count(scene,t,a)/num_trials];

        end

    end

end


%% Results Table

column_names = {
    'Target2Amplitude', 'ROIThreshold_dB', 'CFARAlpha', 'Trials', 'DenseTarget1Hit_pct', 'DenseTarget2Hit_pct', 'SparseTarget1Hit_pct', 'SparseTarget2Hit_pct', 'ROITarget2Coverage_pct', 'DenseConditionedT2Loss_pct', 'MeanROIBins', 'MeanExpandedBins', 'MeanDopplerFFTs', 'MeanCFARCUTs', 'FFTWorkloadReduction_pct', 'DenseOffTargetCellsPerTrial', 'SparseOffTargetCellsPerTrial'
};

results_table = array2table(results_data, "VariableNames", column_names);


%% Output Folder

output_dir = "output/phase2/07_calibrated_cfar_detection_validation";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Save CSV

csv_file = fullfile(output_dir, "original_vs_calibrated_cfar_detection.csv");

writetable(results_table, csv_file);


%% Plot

fig = figure("Position", [100 100 1100 750]);

subplot(2,1,1);

hold on;

for scene = 2:num_scenes

    original_hits = squeeze(sparse_hit_count(scene,:,1,2)) * 100 / num_trials;

    calibrated_hits = squeeze(sparse_hit_count(scene,:,2,2)) * 100 / num_trials;

    plot(threshold_values, original_hits, "-o", "DisplayName", sprintf("Original alpha, amp %.3f", target2_amplitudes(scene)));

    plot(threshold_values, calibrated_hits, "--s", "DisplayName", sprintf("Calibrated alpha, amp %.3f", target2_amplitudes(scene)));

end

xlabel("ROI Threshold [dB]");
ylabel("Target 2 Hit Rate [%]");
title("Weak Target Detection - Original vs Calibrated CFAR");

ylim([0 100]);
grid on;
legend("Location","best");


subplot(2,1,2);

hold on;

for scene = 1:num_scenes

    reduction = 100 * (1 - fft_count_sum(scene,:) / (num_trials*fair_dense_ffts));

    plot(threshold_values, reduction, "-o", "DisplayName", sprintf("Target 2 amp %.3f", target2_amplitudes(scene)));

end

xlabel("ROI Threshold [dB]");
ylabel("Doppler FFT Workload Reduction [%]");
title("Sparse FFT Workload Including CFAR Margin");

ylim([0 100]);
grid on;
legend("Location","best");


plot_file = fullfile(output_dir, "calibrated_cfar_weak_target_tradeoff.png");

exportgraphics(fig, plot_file, "Resolution", 300);


%% Save TXT Summary

txt_file = fullfile(output_dir, "calibrated_cfar_detection_validation.txt");

file = fopen(txt_file, "w");

assert(file ~= -1, "Cannot create Group 3D result file");

fprintf(file, "=== Group 3D: Calibrated CFAR Detection Validation ===\n\n");

fprintf(file, "Trials per amplitude = %d\n", num_trials);
fprintf(file, "Noise SNR = %.1f dB\n", snr_db);
fprintf(file, "Original alpha = %.6f\n", alpha_original);
fprintf(file, "Calibrated alpha = %.6f\n", alpha_calibrated);
fprintf(file, "ROI thresholds = 3 and 6 dB\n");
fprintf(file, "Target matching = +/-1 range bin and +/-1 Doppler bin\n\n");

for scene = 1:num_scenes

    fprintf(file, "[Target 2 amplitude = %.3f]\n", target2_amplitudes(scene));

    for t = 1:num_thresholds

        fprintf(file, "\nROI threshold = %.1f dB\n", threshold_values(t));

        fprintf(file, "Mean ROI bins = %.2f\n", roi_bin_sum(scene,t)/num_trials);
        fprintf(file, "Mean expanded bins = %.2f\n", expanded_bin_sum(scene,t)/num_trials);
        fprintf(file, "Mean Sparse Doppler FFTs = %.2f\n", fft_count_sum(scene,t)/num_trials);

        for a = 1:num_alpha

            fprintf(file, "Alpha = %.6f\n", alpha_values(a));
            fprintf(file, "Dense Target 2 hit rate = %.2f%%\n", 100*dense_hit_count(scene,a,2)/num_trials);
            fprintf(file, "Sparse Target 2 hit rate = %.2f%%\n", 100*sparse_hit_count(scene,t,a,2)/num_trials);
            fprintf(file, "Off-target Dense cells per trial = %.2f\n", dense_offtarget_count(scene,a)/num_trials);
            fprintf(file, "Off-target Sparse cells per trial = %.2f\n", sparse_offtarget_count(scene,t,a)/num_trials);

        end

    end

    fprintf(file, "\n");

end

fprintf(file, "[Validation]\n");
fprintf(file, "Selective FFT and CFAR spot-check = PASS\n");
fprintf(file, "Sparse detection subset of Dense = PASS\n");
fprintf(file, "Calibrated alpha detection monotonicity = PASS\n");
fprintf(file, "Workload counter validation = PASS\n");
fprintf(file, "Overall code validation = PASS\n\n");

fprintf(file, "[Limitations]\n");
fprintf(file, "Sparse detections use Dense CFAR decisions restricted to ROI CUTs.\n");
fprintf(file, "Actual Sparse FFT and CFAR are spot-checked on one trial.\n");
fprintf(file, "Target hits use local ground-truth matching windows.\n");
fprintf(file, "Off-target cells may contain target sidelobes and are not pure-noise false alarms.\n");
fprintf(file, "Noise-only Pfa was validated separately in Group 3C-2.\n");
fprintf(file, "The original CFAR implementation is not modified.\n");
fprintf(file, "FFT reduction is a workload proxy, not measured latency or power.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 3D Validation ---\n");

fprintf("Trials per amplitude = %d\n", num_trials);

fprintf("Original alpha = %.6f\n", alpha_original);

fprintf("Calibrated alpha = %.6f\n", alpha_calibrated);

fprintf("Selective FFT and CFAR regression = PASS\n");

fprintf("Sparse detection subset = PASS\n");

fprintf("Calibrated alpha detection monotonicity = PASS\n");

fprintf("Group 3D code validation passed.\n");

fprintf("CSV saved: %s\n", csv_file);
fprintf("Plot saved: %s\n", plot_file);
fprintf("TXT saved: %s\n", txt_file);
