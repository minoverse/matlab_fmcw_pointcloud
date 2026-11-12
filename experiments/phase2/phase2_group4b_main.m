function phase2_group4b_main()
%% Phase 2 Group 4B - Bounded Workload Control

rng(1);

base = phase2_main();
cfg = base.cfg;
c = base.c;
fc = base.fc;
B = base.B;
Tchirp = base.Tchirp;
fs = base.fs;
Nchirp = base.Nchirp;
Ns = base.Ns;
half_N = base.half_N;
S = base.S;
lambda = base.lambda;
Nrx = base.Nrx;
Nrx_az = base.Nrx_az;
Nrx_el = base.Nrx_el;
d_az = base.d_az;
d_el = base.d_el;
rx_y = base.rx_y;
rx_z = base.rx_z;
NangleFFT_az = base.NangleFFT_az;
NangleFFT_el = base.NangleFFT_el;
target_range = base.target_range;
target_velocity = base.target_velocity;
target_angle = base.target_angle;
target_elevation = base.target_elevation;
target_amplitude = base.target_amplitude;
Ntargets = base.Ntargets;
range_resolution = base.range_resolution;
tau = base.tau;
fb = base.fb;
fd = base.fd;
Tr = base.Tr;
Td = base.Td;
Gr = base.Gr;
Gd = base.Gd;
Pfa = base.Pfa;
snr_db = base.snr_db;
adc_cube = base.adc_cube;
t_fast = base.t_fast;
range_data = base.range_data;
range_axis = base.range_axis;
range_data_no_window = base.range_data_no_window;
roi_threshold_db = base.roi_threshold_db;
roi = base.roi;
range_positive = base.range_positive;
range_energy = base.range_energy;
roi_noise_floor = base.roi_noise_floor;
roi_threshold = base.roi_threshold;
roi_mask = base.roi_mask;
roi_bins = base.roi_bins;
roi_count = base.roi_count;
roi_fraction = base.roi_fraction;
target_roi_selected = base.target_roi_selected;
apply_clutter_removal = base.apply_clutter_removal;
range_data_for_doppler = base.range_data_for_doppler;
range_doppler_rx = base.range_doppler_rx;
RD_power = base.RD_power;
velocity_axis = base.velocity_axis;
RD_mag = base.RD_mag;
RD_db = base.RD_db;
cfar_map = base.cfar_map;
noise_map = base.noise_map;
det_range_idx = base.det_range_idx;
det_doppler_idx = base.det_doppler_idx;
cfar_stats = base.cfar_stats;
Ntraining = base.Ntraining;
alpha = base.alpha;
total_cells = base.total_cells;
guard_cells = base.guard_cells;
peak_range_idx = base.peak_range_idx;
peak_doppler_idx = base.peak_doppler_idx;
peak_power = base.peak_power;
Npeaks = base.Npeaks;
estimated_azimuth = base.estimated_azimuth;
estimated_elevation = base.estimated_elevation;
estimated_snr = base.estimated_snr;
data = base.data;
cost = base.cost;
range_fft_count = base.range_fft_count;
range_fft_diagnostic_count = base.range_fft_diagnostic_count;
doppler_fft_current = base.doppler_fft_current;
doppler_fft_fair = base.doppler_fft_fair;
doppler_fft_sparse = base.doppler_fft_sparse;
projected_roi_reduction = base.projected_roi_reduction;
cfar_cut_count = base.cfar_cut_count;
cfar_training_accesses = base.cfar_training_accesses;


addpath("tests");

assert(isfile("src/optimization/roi_merge_blocks.m"));
assert(isfile("src/optimization/block_doppler_processing.m"));
assert(isfile("src/optimization/roi_budget_selection.m"));


%% Settings

scene_amplitudes = [0.7, 0.005];

roi_thresholds = [3, 6];

profile_names = ["Strict", "Balanced", "Relaxed", "Unbounded"];

top_k_limits = [3, 6, 12, half_N];

fft_caps = [100, 200, 400, half_N*Nrx];

num_scenes = length(scene_amplitudes);
num_thresholds = length(roi_thresholds);
num_profiles = length(profile_names);

margin = Tr + Gr;

fair_dense_ffts = half_N * Nrx;

assert(fair_dense_ffts == 2000);


%% Range and Doppler Windows

range_n = (0:Ns-1)';

range_window_b = 0.5 - 0.5*cos(2*pi*range_n/(Ns-1));

doppler_n = (0:Nchirp-1)';

doppler_window_b = 0.5 - 0.5*cos(2*pi*doppler_n/Nchirp);


%% Dense CFAR Kernel For Weak Scene

doppler_margin = Td + Gd;

cfar_kernel = ones(2*margin+1, 2*doppler_margin+1);

cfar_kernel(margin+1-Gr:margin+1+Gr, doppler_margin+1-Gd:doppler_margin+1+Gd) = 0;

assert(nnz(cfar_kernel) == Ntraining);

valid_range = margin+1:half_N-margin;

valid_doppler = doppler_margin+1:Nchirp-doppler_margin;


%% Results

num_rows = num_scenes * num_thresholds * num_profiles;

results_data = zeros(num_rows,17);

profile_labels = strings(num_rows,1);

result_row = 0;

regression_count = 0;


%% Scene Loop

for scene = 1:num_scenes

    target_amplitude(2) = scene_amplitudes(scene);

    fprintf("\n--- Group 4B Scene: Target 2 amplitude %.3f ---\n", target_amplitude(2));


    %% Scene 1 Uses Existing Dense Baseline

    if scene == 2

        %% New Weak-Target Scene

        cfg_scene = radar_config();
cfg_scene.target_amplitude = target_amplitude;
[adc_cube, t_fast, noise] = simulate_adc_cube_process(cfg_scene, 20);

        % ADC Simulation process does not create a figure

        adc_windowed_b = adc_cube .* reshape(range_window_b, Ns, 1, 1);

        range_data = fft(adc_windowed_b, [], 1);

        range_data_for_doppler = range_data;


        %% Fair Dense Positive-Range Doppler FFT

        dense_input = range_data_for_doppler(1:half_N,:,:) .* reshape(doppler_window_b,1,1,Nchirp);

        range_doppler_rx = fft(dense_input,Nchirp,3);

        range_doppler_rx = fftshift(range_doppler_rx,3);


        %% Optional Static Clutter Notch

        if apply_clutter_removal

            center_bin = floor(Nchirp/2)+1;

            range_doppler_rx(:,:,center_bin-1:center_bin+1) = 0;

        end


        %% Dense RD Power

        RD_power = reshape(sum(abs(range_doppler_rx).^2,2),half_N,Nchirp);


        %% Dense CA-CFAR

        training_sum = conv2(RD_power,cfar_kernel,"same");

        noise_map = zeros(half_N,Nchirp);

        noise_map(valid_range,valid_doppler) = training_sum(valid_range,valid_doppler)/Ntraining;

        cfar_map = false(half_N,Nchirp);

        cfar_map(valid_range,valid_doppler) = RD_power(valid_range,valid_doppler) > alpha*noise_map(valid_range,valid_doppler);

    end


    %% Positive Range FFT and Energy

    range_positive = range_data(1:half_N,:,:);

    range_energy_b = squeeze(sum(sum(abs(range_positive).^2,3),2));

    noise_floor_b = median(range_energy_b);

    dense_fft = range_doppler_rx(1:half_N,:,:);


    %% Dense Peak Grouping Reference

    [dense_det_r, dense_det_d] = find(cfar_map);

    dense_peak_pairs = zeros(0,2);

    for p = 1:length(dense_det_r)

        r = dense_det_r(p);
        d = dense_det_d(p);

        if r <= 1 || r >= half_N || d <= 1 || d >= Nchirp
            continue;
        end

        local_area = RD_power(r-1:r+1,d-1:d+1);

        if RD_power(r,d) == max(local_area(:))

            dense_peak_pairs(end+1,:) = [r,d];

        end

    end


    %% Ground Truth Bins

    target_r = zeros(Ntargets,1);
    target_d = zeros(Ntargets,1);

    for k = 1:Ntargets

        [~,target_r(k)] = min(abs(range_axis-target_range(k)));

        [~,target_d(k)] = min(abs(velocity_axis-target_velocity(k)));

    end


    %% ROI Threshold Loop

    for t = 1:num_thresholds

        threshold_db = roi_thresholds(t);

        threshold_value = noise_floor_b * 10^(threshold_db/10);

        candidate_bins = find(range_energy_b > threshold_value);


        %% Budget Profiles

        for p = 1:num_profiles

            max_rois = top_k_limits(p);

            max_fft_count = fft_caps(p);


            %% Bounded Selection

            [selected_bins, roi_blocks, expanded_bins, dropped_bins] = roi_budget_selection(candidate_bins,range_energy_b,half_N,margin,max_rois,max_fft_count,Nrx);


            %% Actual Block Doppler Processing

            [block_doppler_rx, block_RD_power, block_fft_count] = block_doppler_processing(range_data_for_doppler,roi_blocks,apply_clutter_removal);


            %% Hard Workload Validation

            assert(length(selected_bins) <= max_rois);

            assert(block_fft_count <= max_fft_count);

            assert(block_fft_count == length(expanded_bins)*Nrx);

            assert(length(selected_bins)+length(dropped_bins) == length(candidate_bins));


            %% Dense FFT Regression

            if ~isempty(expanded_bins)

                dense_selected_fft = dense_fft(expanded_bins,:,:);

                fft_error = max(abs(block_doppler_rx(:)-dense_selected_fft(:)));

                fft_tolerance = 1e-10*max(1,max(abs(dense_selected_fft(:))));

                assert(fft_error <= fft_tolerance);

                dense_selected_power = RD_power(expanded_bins,:);

                power_error = max(abs(block_RD_power(:)-dense_selected_power(:)));

                power_tolerance = 1e-10*max(1,max(abs(dense_selected_power(:))));

                assert(power_error <= power_tolerance);

            end


            %% ROI-Only CFAR

            roi_bins = selected_bins;

            active_range_bins = expanded_bins;

            sparse_RD_power = block_RD_power;

            radar_roi = struct("half_N",half_N,"Nchirp",Nchirp,"Tr",Tr,"Td",Td,"Gr",Gr,"Gd",Gd,"Pfa",Pfa,"Ntraining",Ntraining);
[roi_cfar_local, roi_noise_local, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, roi_stats, roi_det_range_idx, roi_det_doppler_idx, range_lookup] = roi_only_cfar_process(roi_bins, active_range_bins, sparse_RD_power, radar_roi);
print_roi_only_cfar(roi_stats);


            %% Dense vs Bounded CFAR Regression

            expected_cfar_map = cfar_map(roi_bins,:);

            assert(isequal(roi_cfar_local,expected_cfar_map));

            if ~isempty(roi_bins)

                dense_selected_noise = noise_map(roi_bins,:);

                noise_error = max(abs(roi_noise_local(:)-dense_selected_noise(:)));

                noise_tolerance = 1e-10*max(1,max(abs(dense_selected_noise(:))));

                assert(noise_error <= noise_tolerance);

            end


            %% Dense vs Bounded Peak Grouping

            expected_peaks = dense_peak_pairs(ismember(dense_peak_pairs(:,1),roi_bins),:);

            bounded_peaks = [roi_peak_range_idx(:),roi_peak_doppler_idx(:)];

            assert(isequal(sortrows(expected_peaks),sortrows(bounded_peaks)));


            %% Compare Target Hits

            bounded_cfar_map = false(half_N,Nchirp);

            bounded_cfar_map(roi_bins,:) = roi_cfar_local;

            dense_hit = zeros(Ntargets,1);

            sparse_hit = zeros(Ntargets,1);

            for k = 1:Ntargets

                r0 = target_r(k);
                d0 = target_d(k);

                r_near = max(1,r0-1):min(half_N,r0+1);

                d_near = max(1,d0-1):min(Nchirp,d0+1);

                dense_hit(k) = any(any(cfar_map(r_near,d_near)));

                sparse_hit(k) = any(any(bounded_cfar_map(r_near,d_near)));

                assert(sparse_hit(k) <= dense_hit(k));

            end


            %% Save Results

            result_row = result_row+1;

            profile_labels(result_row) = profile_names(p);

            fft_reduction = 100*(1-block_fft_count/fair_dense_ffts);

            results_data(result_row,:) = [scene_amplitudes(scene),threshold_db,max_rois,max_fft_count,length(candidate_bins),length(selected_bins),length(dropped_bins),size(roi_blocks,1),length(expanded_bins),block_fft_count,fft_reduction,roi_cfar_cut_count,nnz(roi_cfar_local),dense_hit(1),dense_hit(2),sparse_hit(1),sparse_hit(2)];

            regression_count = regression_count+1;


            %% Print

            fprintf("%s | ROI %.1f dB | Candidates %d | Selected %d | FFT %d/%d | Target 2 %d/%d\n",profile_names(p),threshold_db,length(candidate_bins),length(selected_bins),block_fft_count,max_fft_count,sparse_hit(2),dense_hit(2));

        end

    end

end


%% Synthetic Stress Tests

[stress_rows,stress_candidate_counts,stress_fft_caps] = roi_budget_regression(half_N,Nrx,Tr,Gr);


%% Prepare Results Table

column_names = {'Target2Amplitude','ROIThreshold_dB','TopKLimit','FFTCap','CandidateROIBins','SelectedROIBins','DroppedROIBins','ProcessingBlocks','ExpandedRangeBins','ActualDopplerFFTs','FFTWorkloadReduction_pct','CFARCUTs','CFARDetectionCells','DenseTarget1Hit','DenseTarget2Hit','BoundedTarget1Hit','BoundedTarget2Hit'};

results_table = array2table(results_data,"VariableNames",column_names);

results_table.Profile = profile_labels;


%% Output Folder

output_dir = "output/phase2/09_topk_bounded_workload";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Save Real Scene CSV

results_csv = fullfile(output_dir,"topk_bounded_detection_results.csv");

writetable(results_table,results_csv);


%% Save Synthetic Stress CSV

stress_column_names = {'CandidateROIBins','FFTCap','SelectedROIBins','ProcessingBlocks','ExpandedRangeBins','ProjectedDopplerFFTs'};

stress_table = array2table(stress_rows,"VariableNames",stress_column_names);

stress_csv = fullfile(output_dir,"bounded_workload_density_stress.csv");

writetable(stress_table,stress_csv);


%% Plot

% 이전 Figure 정리
close all;

fig = figure("Position",[100 100 1100 1050]);


%% Weak Target FFT Budget

weak_rows = find(abs(results_data(:,1)-0.005)<1e-12 & results_data(:,2)==3);

assert(length(weak_rows)==num_profiles);

subplot(3,1,1);

plot(1:num_profiles,results_data(weak_rows,10),"-o","LineWidth",1.5);

hold on;

plot(1:num_profiles,results_data(weak_rows,4),"--s","LineWidth",1.2);

xticks(1:num_profiles);
xticklabels(profile_names);

xlabel("Workload Profile");
ylabel("Doppler FFT Count");

title("Weak Target Scene: Actual FFT Workload vs Budget");

legend("Actual FFTs","FFT Hard Cap","Location","best");

grid on;


%% Weak Target Detection

subplot(3,1,2);

plot(1:num_profiles,results_data(weak_rows,17),"-o","LineWidth",1.5);

hold on;

plot(1:num_profiles,results_data(weak_rows,15),"--s","LineWidth",1.2);

xticks(1:num_profiles);
xticklabels(profile_names);

xlabel("Workload Profile");
ylabel("Target 2 Hit");

title("Weak Target Detection: Bounded vs Dense CFAR");

legend("Bounded ROI CFAR","Dense CFAR","Location","best");

ylim([-0.1 1.1]);

grid on;


%% Synthetic Density Stress

subplot(3,1,3);

hold on;

plot(stress_rows(stress_rows(:,2)==100,1),stress_rows(stress_rows(:,2)==100,6),"-o","LineWidth",1.3);

plot(stress_rows(stress_rows(:,2)==200,1),stress_rows(stress_rows(:,2)==200,6),"-s","LineWidth",1.3);

yline(100,"--");
yline(200,"--");

xlabel("Candidate ROI Bin Count");
ylabel("Projected Doppler FFTs");

title("Synthetic ROI Density Stress: Hard Workload Caps");

legend("Cap 100","Cap 200","100 FFT Limit","200 FFT Limit","Location","best");

grid on;


%% Save PNG
plot_file = fullfile(output_dir, "topk_budget_and_weak_target_tradeoff.png");

plot_saved = false;
try
    drawnow;
    assert(isgraphics(fig, "figure"), "Figure is invalid or closed");
    frame = getframe(fig);
    assert(~isempty(frame.cdata), "Captured frame is empty");
    imwrite(frame.cdata, plot_file);
    plot_saved = isfile(plot_file);
catch ME
    warning("PNG capture failed: %s", ME.message);
end

%% Save TXT

txt_file = fullfile(output_dir,"topk_bounded_workload_validation.txt");

file = fopen(txt_file,"w");

assert(file ~= -1,"Cannot create Group 4B summary");

fprintf(file,"=== Group 4B: Top-K and Bounded Workload Control ===\n\n");

fprintf(file,"[Settings]\n");
fprintf(file,"CFAR range margin = %d\n",margin);
fprintf(file,"Fair Dense Doppler FFTs = %d\n",fair_dense_ffts);
fprintf(file,"CFAR alpha = %.6f\n",alpha);
fprintf(file,"Scene amplitudes = 0.7 and 0.005\n");
fprintf(file,"ROI thresholds = 3 and 6 dB\n");
fprintf(file,"Profiles = Strict, Balanced, Relaxed, Unbounded\n\n");

fprintf(file,"[Real Scene Tests]\n");

for i = 1:result_row

    fprintf(file,"Amplitude %.3f | ROI %.1f dB | %s | Candidates %d | Selected %d | FFTs %d/%d | Dense T2 %d | Bounded T2 %d\n",results_data(i,1),results_data(i,2),char(profile_labels(i)),results_data(i,5),results_data(i,6),results_data(i,10),results_data(i,4),results_data(i,15),results_data(i,17));

end

fprintf(file,"\n[Synthetic Stress]\n");
fprintf(file,"Candidate densities = %d\n",length(stress_candidate_counts));
fprintf(file,"FFT cap settings = %d\n",length(stress_fft_caps));
fprintf(file,"Total stress cases = %d\n",size(stress_rows,1));

fprintf(file,"\n[Validation]\n");
fprintf(file,"Top-K hard limit = PASS\n");
fprintf(file,"Actual Doppler FFT hard limit = PASS\n");
fprintf(file,"CFAR margin coverage = PASS\n");
fprintf(file,"Dense vs Block FFT regression = PASS\n");
fprintf(file,"Dense vs Bounded CFAR regression = PASS\n");
fprintf(file,"Dense vs Bounded peak grouping = PASS\n");
fprintf(file,"Empty and full budget tests = PASS\n");
fprintf(file,"Overall code validation = PASS\n");

fprintf(file,"\n[Limitations]\n");
fprintf(file,"ROI ranking uses noncoherent range energy only.\n");
fprintf(file,"Greedy ranking is not guaranteed to preserve weak targets.\n");
fprintf(file,"Weak target test uses one noise realization per scene.\n");
fprintf(file,"Synthetic stress tests verify selection counters, not measured FFT execution.\n");
fprintf(file,"Doppler FFT budget excludes Range FFT, CFAR and AoA workload.\n");
fprintf(file,"FFT reduction is a workload proxy, not measured latency or power.\n");
fprintf(file,"The original CFAR implementation is unchanged.\n");

fclose(file);


%% Final Validation

assert(result_row==num_rows);
assert(regression_count==num_rows);
assert(all(results_data(:,10)<=results_data(:,4)));
assert(all(results_data(:,6)<=results_data(:,3)));


fprintf("\n--- Group 4B Validation ---\n");

fprintf("Real scene configurations tested = %d\n",result_row);

fprintf("Synthetic stress configurations tested = %d\n",size(stress_rows,1));

fprintf("Top-K hard limit = PASS\n");

fprintf("Actual Doppler FFT hard limit = PASS\n");

fprintf("Dense FFT / CFAR / Peak regression = PASS\n");

fprintf("Synthetic ROI density stress = PASS\n");

fprintf("Group 4B code validation passed.\n");

fprintf("Results CSV: %s\n",results_csv);

fprintf("Stress CSV: %s\n",stress_csv);

if plot_saved
    fprintf("Plot: %s\n", plot_file);
else
    fprintf("Plot: NOT SAVED (graphics issue)\n");
end

fprintf("Summary: %s\n",txt_file);
end
