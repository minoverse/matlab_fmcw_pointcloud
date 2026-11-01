
%% Phase 2 Group 5C - Final Statistical Validation

clear;
clc;

assert(isfile("config.m"),"Run from the project root");

addpath(genpath("src"));

run("config.m");

assert(Ntargets==2,"This experiment expects two targets");


%% Experiment Settings

snr_values = [0,10,20,30];

weak_amplitudes = [0.03,0.005];

roi_thresholds = [2,3,4,5,6];

alpha_values = [216*(1e-5^(-1/216)-1),4.8];

method_names = ["Dense","Sparse ROI","Top-K Bins","Top-K Blocks"];

num_trials = 200;

max_fft_count = 200;

top_k_limit = 6;


%% CFAR Settings

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

margin = Tr+Gr;
doppler_margin = Td+Gd;

half_N = floor(Ns/2);

Ntraining = (2*margin+1)*(2*doppler_margin+1)-(2*Gr+1)*(2*Gd+1);

assert(num_trials>=1);
assert(Ntraining==216);
assert(half_N*Nrx==2000);


%% Result Dimensions

num_scenes = length(snr_values)*length(weak_amplitudes);

num_thresholds = length(roi_thresholds);

num_methods = length(method_names);

num_alpha = length(alpha_values);


%% Result Counters

hit_counts = zeros(num_scenes,num_thresholds,num_methods,num_alpha,2);

coverage_counts = zeros(num_scenes,num_thresholds,num_methods,2);

fft_sum = zeros(num_scenes,num_thresholds,num_methods);

fft_max = zeros(num_scenes,num_thresholds,num_methods);

cfar_cut_sum = zeros(num_scenes,num_thresholds,num_methods);

scene_complete = false(num_scenes,1);


%% Output Folder

output_dir = "output/phase2/14_final_statistical_validation";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Checkpoint For Interrupted Runs

checkpoint_file = fullfile(output_dir,"group5c_checkpoint.mat");

if isfile(checkpoint_file)

    previous = load(checkpoint_file);

    assert(isequal(previous.snr_values,snr_values),"Checkpoint SNR settings differ");

    assert(isequal(previous.weak_amplitudes,weak_amplitudes),"Checkpoint amplitudes differ");

    assert(isequal(previous.roi_thresholds,roi_thresholds),"Checkpoint thresholds differ");

    assert(isequal(previous.alpha_values,alpha_values),"Checkpoint alpha values differ");

    assert(previous.num_trials==num_trials,"Checkpoint trial count differs");

    hit_counts = previous.hit_counts;

    coverage_counts = previous.coverage_counts;

    fft_sum = previous.fft_sum;

    fft_max = previous.fft_max;

    cfar_cut_sum = previous.cfar_cut_sum;

    scene_complete = previous.scene_complete;

    fprintf("Loaded checkpoint: %d/%d scenes already complete\n",nnz(scene_complete),num_scenes);

end


%% FFT Windows

n_range = (0:Ns-1)';

range_window = 0.5-0.5*cos(2*pi*n_range/(Ns-1));

n_doppler = (0:Nchirp-1)';

doppler_window = 0.5-0.5*cos(2*pi*n_doppler/Nchirp);


%% CFAR Window

valid_range = margin+1:half_N-margin;

valid_doppler = doppler_margin+1:Nchirp-doppler_margin;

cfar_kernel = ones(2*margin+1,2*doppler_margin+1);

cfar_kernel(margin+1-Gr:margin+1+Gr,doppler_margin+1-Gd:doppler_margin+1+Gd) = 0;

assert(nnz(cfar_kernel)==Ntraining);


%% Range and Velocity Axes

range_axis = c*(0:half_N-1)*fs/(2*S*Ns);

velocity_axis = (-Nchirp/2:Nchirp/2-1)*(1/Tchirp)*lambda/(2*Nchirp);


%% Ground Truth Matching Bins

range_rows = cell(2,1);

doppler_cols = cell(2,1);

for k = 1:2

    [~,r0] = min(abs(range_axis-target_range(k)));

    [~,d0] = min(abs(velocity_axis-target_velocity(k)));

    range_rows{k} = max(1,r0-1):min(half_N,r0+1);

    doppler_cols{k} = max(1,d0-1):min(Nchirp,d0+1);

end


%% Monte Carlo

fprintf("\n--- Group 5C Final Monte Carlo ---\n");

fprintf("Scenes = %d, trials per scene = %d, paired thresholds = %d\n",num_scenes,num_trials,num_thresholds);

scene_index = 0;

for s = 1:length(snr_values)

    for a = 1:length(weak_amplitudes)

        scene_index = scene_index+1;


        %% Skip Completed Scene

        if scene_complete(scene_index)

            fprintf("Scene %d/%d already complete - skipping\n",scene_index,num_scenes);

            continue;

        end


        %% Independent Random Seed

        rng(5000+100*s+a);

        snr_db = snr_values(s);

        amplitudes = [1,weak_amplitudes(a)];


        %% Build Clean ADC Signal Once Per Scene

        clean_adc = complex(zeros(Ns,Nrx,Nchirp));

        t_fast = (0:Ns-1)'/fs;

        chirp_numbers = 0:Nchirp-1;

        for k = 1:2

            delay = 2*target_range(k)/c;

            fb_k = S*delay;

            fd_k = 2*target_velocity(k)/lambda;

            fast_wave = exp(1j*2*pi*fb_k*t_fast);

            slow_wave = exp(1j*2*pi*fd_k*chirp_numbers*Tchirp);

            wave = fast_wave*slow_wave;

            direction_y = cosd(target_elevation(k))*sind(target_angle(k));

            direction_z = sind(target_elevation(k));

            for rx = 1:Nrx

                spatial_phase = 2*pi/lambda*(rx_y(rx)*direction_y+rx_z(rx)*direction_z);

                clean_adc(:,rx,:) = clean_adc(:,rx,:) + reshape(amplitudes(k)*exp(1j*spatial_phase)*wave,Ns,1,Nchirp);

            end

        end


        %% Noise Power

        noise_power = max(amplitudes)^2/10^(snr_db/10);


        %% Independent Trials

        for trial = 1:num_trials

            noise = sqrt(noise_power/2)*(randn(Ns,Nrx,Nchirp)+1j*randn(Ns,Nrx,Nchirp));

            adc_cube = clean_adc+noise;


            %% Range FFT

            range_data = fft(adc_cube.*reshape(range_window,Ns,1,1),[],1);

            range_positive = range_data(1:half_N,:,:);

            range_energy = squeeze(sum(sum(abs(range_positive).^2,3),2));

            noise_floor = median(range_energy);


            %% Dense Doppler FFT

            doppler_input = range_positive.*reshape(doppler_window,1,1,Nchirp);

            dense_fft = fftshift(fft(doppler_input,Nchirp,3),3);

            RD_power = reshape(sum(abs(dense_fft).^2,2),half_N,Nchirp);


            %% Faster Separable CFAR Training Sum

            full_training = conv2(conv2(RD_power,ones(2*margin+1,1),"same"),ones(1,2*doppler_margin+1),"same");

            guard_training = conv2(conv2(RD_power,ones(2*Gr+1,1),"same"),ones(1,2*Gd+1),"same");

            training_sum = full_training-guard_training;


            %% Verify Against Original CFAR Kernel

            if trial==1

                reference_sum = conv2(RD_power,cfar_kernel,"same");

                training_error = max(abs(reference_sum(:)-training_sum(:)));

                training_tolerance = 1e-10*max(1,max(abs(reference_sum(:))));

                assert(training_error<=training_tolerance,"Fast CFAR training mismatch");

            end

            noise_map = training_sum/Ntraining;


            %% Dense CFAR For Both Alpha Values

            dense_cfar = false(half_N,Nchirp,num_alpha);

            for q = 1:num_alpha

                map = false(half_N,Nchirp);

                map(valid_range,valid_doppler) = RD_power(valid_range,valid_doppler)>alpha_values(q)*noise_map(valid_range,valid_doppler);

                dense_cfar(:,:,q) = map;

            end


            %% Paired ROI Thresholds

            for t = 1:num_thresholds

                threshold_db = roi_thresholds(t);

                threshold = noise_floor*10^(threshold_db/10);

                candidate_bins = find(range_energy>threshold);


                %% Sparse ROI Without Top-K Limit

                [~,sparse_expanded] = roi_merge_blocks(candidate_bins,half_N,margin);


                %% Top-K Bins

                [bin_selected,~,bin_expanded,~] = roi_budget_selection(candidate_bins,range_energy,half_N,margin,top_k_limit,max_fft_count,Nrx);


                %% Top-K Blocks

                [block_selected,~,~,block_expanded,~,~] = roi_block_budget_selection(candidate_bins,range_energy,half_N,margin,top_k_limit,max_fft_count,Nrx);


                %% Selected Range CUT Masks

                selected_mask = false(half_N,num_methods);

                selected_mask(:,1) = true;

                selected_mask(candidate_bins,2) = true;

                selected_mask(bin_selected,3) = true;

                selected_mask(block_selected,4) = true;


                %% Doppler FFT Workload

                fft_counts = [half_N*Nrx,length(sparse_expanded)*Nrx,length(bin_expanded)*Nrx,length(block_expanded)*Nrx];

                assert(all(fft_counts<=half_N*Nrx));

                assert(all(fft_counts(3:4)<=max_fft_count));


                %% Store Workload Counters

                for m = 1:num_methods

                    count = fft_counts(m);

                    cfar_cuts = nnz(selected_mask(valid_range,m))*length(valid_doppler);

                    fft_sum(scene_index,t,m) = fft_sum(scene_index,t,m)+count;

                    fft_max(scene_index,t,m) = max(fft_max(scene_index,t,m),count);

                    cfar_cut_sum(scene_index,t,m) = cfar_cut_sum(scene_index,t,m)+cfar_cuts;


                    %% Target ROI Coverage

                    for k = 1:2

                        r_near = range_rows{k};

                        coverage_counts(scene_index,t,m,k) = coverage_counts(scene_index,t,m,k)+double(any(selected_mask(r_near,m)));

                    end

                end


                %% Detection Using Identical Dense CFAR Decisions

                for q = 1:num_alpha

                    map = dense_cfar(:,:,q);

                    for k = 1:2

                        r_near = range_rows{k};

                        d_near = doppler_cols{k};

                        detected_rows = any(map(r_near,d_near),2);

                        dense_hit = any(detected_rows);

                        for m = 1:num_methods

                            hit = any(detected_rows & selected_mask(r_near,m));

                            assert(~hit || dense_hit,"Sparse hit cannot exceed Dense");

                            hit_counts(scene_index,t,m,q,k) = hit_counts(scene_index,t,m,q,k)+double(hit);

                        end

                    end

                end

            end

        end


        %% Save Completed Scene Checkpoint

        scene_complete(scene_index) = true;

        save(checkpoint_file,"snr_values","weak_amplitudes","roi_thresholds","alpha_values","num_trials","hit_counts","coverage_counts","fft_sum","fft_max","cfar_cut_sum","scene_complete");

        fprintf("Scene %d/%d complete: SNR %g dB, target2 amp %.3f, %d trials\n",scene_index,num_scenes,snr_db,amplitudes(2),num_trials);

    end

end

assert(all(scene_complete),"Some scenes not completed");


%% Statistical Results

num_rows = num_scenes*num_thresholds*num_methods*num_alpha;

SNR_dB = zeros(num_rows,1);
Target2Amplitude = zeros(num_rows,1);
ROIThreshold_dB = zeros(num_rows,1);
CFARAlpha = zeros(num_rows,1);
Method = strings(num_rows,1);

Trials = num_trials*ones(num_rows,1);

Target1Hit_pct = zeros(num_rows,1);
Target2Hit_pct = zeros(num_rows,1);
Target2Miss_pct = zeros(num_rows,1);

Target2HitCI95Lower_pct = zeros(num_rows,1);
Target2HitCI95Upper_pct = zeros(num_rows,1);

Target2ROICoverage_pct = zeros(num_rows,1);

DenseConditionedTarget2Loss_pct = NaN(num_rows,1);

MeanDopplerFFTs = zeros(num_rows,1);
MaxDopplerFFTs = zeros(num_rows,1);
MeanCFARCUTs = zeros(num_rows,1);

DopplerReduction_pct = zeros(num_rows,1);
TotalFFTReduction_pct = zeros(num_rows,1);


%% Common FFT Work Proxy

range_fft_work = Nrx*Nchirp*Ns*log2(Ns);

work_per_doppler = Nchirp*log2(Nchirp);

dense_total_fft_work = range_fft_work+half_N*Nrx*work_per_doppler;


%% Build Result Rows

row = 0;

for s = 1:length(snr_values)

    for a = 1:length(weak_amplitudes)

        scene = (s-1)*length(weak_amplitudes)+a;

        for t = 1:num_thresholds

            for m = 1:num_methods

                for q = 1:num_alpha

                    row = row+1;

                    count1 = hit_counts(scene,t,m,q,1);
                    count2 = hit_counts(scene,t,m,q,2);

                    dense_count2 = hit_counts(scene,t,1,q,2);

                    assert(count1<=hit_counts(scene,t,1,q,1));
                    assert(count2<=dense_count2);


                    %% Hit Rate and 95% Wilson Confidence Interval

                    p = count2/num_trials;

                    z = 1.96;

                    denominator = 1+z^2/num_trials;

                    center = (p+z^2/(2*num_trials))/denominator;

                    interval = z*sqrt(p*(1-p)/num_trials+z^2/(4*num_trials^2))/denominator;


                    %% Mean Workload

                    mean_fft = fft_sum(scene,t,m)/num_trials;

                    mean_cuts = cfar_cut_sum(scene,t,m)/num_trials;


                    %% Store Results

                    SNR_dB(row) = snr_values(s);

                    Target2Amplitude(row) = weak_amplitudes(a);

                    ROIThreshold_dB(row) = roi_thresholds(t);

                    CFARAlpha(row) = alpha_values(q);

                    Method(row) = method_names(m);

                    Target1Hit_pct(row) = 100*count1/num_trials;

                    Target2Hit_pct(row) = 100*p;

                    Target2Miss_pct(row) = 100*(1-p);

                    Target2HitCI95Lower_pct(row) = 100*max(0,center-interval);

                    Target2HitCI95Upper_pct(row) = 100*min(1,center+interval);

                    Target2ROICoverage_pct(row) = 100*coverage_counts(scene,t,m,2)/num_trials;

                    MeanDopplerFFTs(row) = mean_fft;

                    MaxDopplerFFTs(row) = fft_max(scene,t,m);

                    MeanCFARCUTs(row) = mean_cuts;

                    DopplerReduction_pct(row) = 100*(1-mean_fft/(half_N*Nrx));

                    TotalFFTReduction_pct(row) = 100*(1-(range_fft_work+mean_fft*work_per_doppler)/dense_total_fft_work);

                    if dense_count2>0

                        DenseConditionedTarget2Loss_pct(row) = 100*(dense_count2-count2)/dense_count2;

                    end


                    %% Hard Cap Validation

                    if m>=3

                        assert(MaxDopplerFFTs(row)<=max_fft_count);

                    end

                end

            end

        end

    end

end


%% Statistical Validation

assert(row==num_rows);

assert(all(Target2Hit_pct>=0 & Target2Hit_pct<=100));

assert(all(Target2HitCI95Lower_pct<=Target2Hit_pct));

assert(all(Target2HitCI95Upper_pct>=Target2Hit_pct));


%% Save CSV

result_table = table(SNR_dB,Target2Amplitude,ROIThreshold_dB,CFARAlpha,Method,Trials,Target1Hit_pct,Target2Hit_pct,Target2Miss_pct,Target2HitCI95Lower_pct,Target2HitCI95Upper_pct,Target2ROICoverage_pct,DenseConditionedTarget2Loss_pct,MeanDopplerFFTs,MaxDopplerFFTs,MeanCFARCUTs,DopplerReduction_pct,TotalFFTReduction_pct);

results_csv = fullfile(output_dir,"final_monte_carlo_detection_workload.csv");

summary_txt = fullfile(output_dir,"final_statistical_validation.txt");

writetable(result_table,results_csv);


%% TXT Summary

file = fopen(summary_txt,"w");

assert(file~=-1,"Cannot save Group 5C summary");

fprintf(file,"=== Group 5C: Statistical Validation ===\n\n");

fprintf(file,"Independent noise trials per SNR/amplitude scene = %d\n",num_trials);

fprintf(file,"SNR values = 0, 10, 20, 30 dB\n");

fprintf(file,"Target 2 amplitudes = 0.03 and 0.005\n");

fprintf(file,"ROI thresholds = 2, 3, 4, 5, 6 dB\n");

fprintf(file,"CFAR alphas = %.6f and %.6f\n",alpha_values(1),alpha_values(2));

fprintf(file,"Top-K Bins/Blocks limit = %d; FFT cap = %d\n",top_k_limit,max_fft_count);

fprintf(file,"Independent noise frames = %d\n",num_scenes*num_trials);

fprintf(file,"Statistical result rows = %d\n",num_rows);

fprintf(file,"95%% Wilson confidence intervals are computed per local-CFAR target hit rate.\n");


%% Weak Target Summary

fprintf(file,"\n[Weak Target 0.005, SNR 20 dB, ROI 3 dB]\n");

focus = find(SNR_dB==20 & abs(Target2Amplitude-0.005)<1e-12 & ROIThreshold_dB==3);

for i = 1:length(focus)

    k = focus(i);

    fprintf(file,"%s | Alpha %.3f | Hit %.1f%% | 95%% CI [%.1f, %.1f] | Mean FFT %.1f | Loss vs Dense %.1f%%\n",Method(k),CFARAlpha(k),Target2Hit_pct(k),Target2HitCI95Lower_pct(k),Target2HitCI95Upper_pct(k),MeanDopplerFFTs(k),DenseConditionedTarget2Loss_pct(k));

end


%% Validation Summary

fprintf(file,"\n[Validation]\n");

fprintf(file,"Independent Monte Carlo noise frames completed = PASS\n");

fprintf(file,"CFAR separable-vs-original kernel = PASS\n");

fprintf(file,"Paired threshold comparison = PASS\n");

fprintf(file,"Sparse target hits subset of Dense = PASS\n");

fprintf(file,"Top-K Bins/Blocks hard FFT cap = PASS\n");

fprintf(file,"Statistical interval computation = PASS\n");

if num_trials>=200

    fprintf(file,"200+ trials per SNR/amplitude scene = PASS\n");

else

    fprintf(file,"200+ trials per SNR/amplitude scene = NOT YET (preliminary run)\n");

end

fprintf(file,"Overall code validation = PASS\n");


%% Important Limitations

fprintf(file,"\n[Limitations]\n");

fprintf(file,"Target hits are local CFAR detections within +/-1 range and Doppler bins.\n");

fprintf(file,"Hit rate can include accidental noise detections near ground truth.\n");

fprintf(file,"Selection methods restrict the same Dense CFAR decisions to selected CUTs.\n");

fprintf(file,"Sparse FFT/CFAR are not re-executed each trial; equivalence was tested in Groups 4B-2/4C.\n");

fprintf(file,"Counts are modeled from selected processing bins; no actual latency or energy measured.\n");

fprintf(file,"The nominal CFAR alpha and noise-only calibrated alpha are both tested.\n");

fprintf(file,"Only the specified synthetic two-target geometry and Gaussian noise are evaluated.\n");

fprintf(file,"Confidence intervals are binomial intervals under independent noise-trial assumptions.\n");

fprintf(file,"No PNG is generated. Completed scenes are saved to a MAT checkpoint.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 5C Validation ---\n");

fprintf("Noise scenes = %d\n",num_scenes);

fprintf("Trials per scene = %d\n",num_trials);

fprintf("Total independent noise frames = %d\n",num_scenes*num_trials);

fprintf("Statistical result rows = %d\n",num_rows);

fprintf("Hard FFT cap / Dense subset / CFAR kernel = PASS\n");

if num_trials>=200

    fprintf("200+ Monte Carlo trials per scene = PASS\n");

else

    fprintf("200+ Monte Carlo trials per scene = PRELIMINARY\n");

end

fprintf("Group 5C code validation passed.\n");

fprintf("Results CSV: %s\n",results_csv);

fprintf("Summary: %s\n",summary_txt);

fprintf("Checkpoint: %s\n",checkpoint_file);
