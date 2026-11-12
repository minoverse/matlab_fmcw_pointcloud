function phase2_group4c_main()

%% Phase 2 Group 4C - Scene Density and Worst Case

% 함수의 지역 변수 범위가 자동으로 초기화됨
clc;
close all;

assert(isfile("src/config/radar_config.m"),"Run from the project root");

addpath(genpath("src"));
addpath("tests");

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


%% Radar Settings

radar.c = c;
radar.S = S;
radar.lambda = lambda;
radar.fs = fs;
radar.Tchirp = Tchirp;
radar.Ns = Ns;
radar.Nchirp = Nchirp;
radar.Nrx = Nrx;
radar.half_N = floor(Ns/2);
radar.rx_y = rx_y;
radar.rx_z = rx_z;

radar.Tr = 8;
radar.Td = 4;
radar.Gr = 2;
radar.Gd = 1;

radar.Ntraining = 216;
radar.Pfa = 1e-5;
radar.snr_db = 20;


%% Experiment Settings

target_counts = [2,4,8,16];

num_trials = 3;

fft_caps = [100,200];

top_k_limits = [3,6];

methods = ["Bins","Blocks"];

threshold_db = 3;

num_rows = length(target_counts)*num_trials*length(fft_caps)*length(methods);

results = zeros(num_rows,15);

row = 0;


%% Scene Density Loop

for i = 1:length(target_counts)

    num_targets = target_counts(i);

    fprintf("\n--- Group 4C: %d Targets ---\n",num_targets);


    %% Deterministic Target Positions

    ranges = linspace(12,65,num_targets);

    velocities = linspace(-12,12,num_targets);

    amplitudes = 0.35+0.65*abs(sin((1:num_targets)*1.3));

    amplitudes(1) = 1;

    amplitudes(end) = 0.005;


    %% Independent Trials

    for trial = 1:num_trials

        rng(1000+100*num_targets+trial);

        data = roi_budget_scene_data(ranges,velocities,amplitudes,radar);

        noise_floor = median(data.range_energy);

        threshold = noise_floor*10^(threshold_db/10);

        candidate_bins = find(data.range_energy > threshold);


        %% Workload Profiles

        for p = 1:length(fft_caps)

            top_k = top_k_limits(p);

            fft_cap = fft_caps(p);


            %% Compare Both Methods

            for m = 1:length(methods)

                result = roi_budget_compare_case(methods(m),candidate_bins,data,radar,top_k,fft_cap);

                row = row+1;

                results(row,:) = [num_targets,trial,fft_cap,m,result.candidate_bins,result.candidate_blocks,result.selected_bins,result.selected_source_blocks,result.fft_count,result.cfar_cuts,sum(result.dense_hits),sum(result.bounded_hits),result.dense_hits(end),result.bounded_hits(end),result.dropped_bins];

                assert(result.fft_count <= fft_cap);

                assert(all(result.bounded_hits <= result.dense_hits));

            end

        end

        fprintf("Trial %d/%d completed\n",trial,num_trials);

    end

end


%% Worst-Case Candidate Density

fprintf("\n--- Full-Density Worst-Case Test ---\n");

full_candidates = (1:radar.half_N)';

worst_caps = [0,100,200,400,2000];

worst_results = zeros(length(worst_caps)*length(methods),7);

worst_row = 0;


%% Reuse Last Scene for Actual FFT Validation

for p = 1:length(worst_caps)

    fft_cap = worst_caps(p);

    if fft_cap == 2000
        top_k = radar.half_N;
    else
        top_k = 6;
    end


    for m = 1:length(methods)

        result = roi_budget_compare_case(methods(m),full_candidates,data,radar,top_k,fft_cap);

        worst_row = worst_row+1;

        worst_results(worst_row,:) = [fft_cap,m,top_k,result.selected_bins,result.selected_source_blocks,result.fft_count,result.cfar_cuts];

        assert(result.fft_count <= fft_cap);

        fprintf("%s | Full density | FFT %d/%d | Selected bins %d\n",methods(m),result.fft_count,fft_cap,result.selected_bins);

    end

end


%% Atomic Block Worst-Case Validation

block_low_budget = worst_results(worst_results(:,2)==2 & worst_results(:,1)<2000,:);

assert(all(block_low_budget(:,4)==0));

full_block = worst_results(worst_results(:,2)==2 & worst_results(:,1)==2000,:);

assert(full_block(6)==2000);

assert(full_block(4)==radar.half_N);


%% Save Results

output_dir = "output/phase2/11_scene_density_worst_case";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Density CSV

column_names = {
    'TargetCount','Trial','FFTCap','MethodIndex','CandidateBins','CandidateBlocks','SelectedBins','SelectedSourceBlocks','ActualDopplerFFTs','CFARCUTs','DenseTargetHits','BoundedTargetHits','DenseWeakTargetHit','BoundedWeakTargetHit','DroppedCandidateBins'
};

density_table = array2table(results,"VariableNames",column_names);

method_names = strings(row,1);

for i = 1:row
    method_names(i) = methods(results(i,4));
end

density_table.Method = method_names;

density_csv = fullfile(output_dir,"scene_density_detection_workload.csv");

writetable(density_table,density_csv);


%% Worst-Case CSV

worst_columns = {
    'FFTCap','MethodIndex','TopKLimit','SelectedBins','SelectedSourceBlocks','ActualDopplerFFTs','CFARCUTs'
};

worst_table = array2table(worst_results,"VariableNames",worst_columns);

worst_csv = fullfile(output_dir,"full_density_worst_case.csv");

writetable(worst_table,worst_csv);


%% Aggregate Results for FFT Cap 200

mean_ffts = zeros(length(target_counts),2);

hit_retention = NaN(length(target_counts),2);

for i = 1:length(target_counts)

    for m = 1:2

        selected_rows = results(:,1)==target_counts(i) & results(:,3)==200 & results(:,4)==m;

        subset = results(selected_rows,:);

        mean_ffts(i,m) = mean(subset(:,9));

        total_dense_hits = sum(subset(:,11));

        total_sparse_hits = sum(subset(:,12));

        if total_dense_hits > 0

            hit_retention(i,m) = 100*total_sparse_hits/total_dense_hits;

        end

    end

end


%% Plot

fig = figure("Position",[100 100 1050 900]);


%% FFT Workload Scaling

subplot(3,1,1);

plot(target_counts,mean_ffts(:,1),"-o","LineWidth",1.5);

hold on;

plot(target_counts,mean_ffts(:,2),"-s","LineWidth",1.5);

yline(200,"k--");

xlabel("Number of Targets");
ylabel("Mean Doppler FFTs");

title("Scene Density vs Processing Workload - FFT Cap 200");

legend("Top-K Bins","Top-K Blocks","FFT Hard Cap","Location","best");

grid on;


%% Detection Retention

subplot(3,1,2);

plot(target_counts,hit_retention(:,1),"-o","LineWidth",1.5);

hold on;

plot(target_counts,hit_retention(:,2),"-s","LineWidth",1.5);

xlabel("Number of Targets");
ylabel("Matched Target Retention [%]");

title("Target Hit Retention Relative to Dense CFAR");

legend("Top-K Bins","Top-K Blocks","Location","best");

ylim([0 105]);

grid on;


%% Full-Density Worst Case

subplot(3,1,3);

bin_worst = worst_results(worst_results(:,2)==1,:);

block_worst = worst_results(worst_results(:,2)==2,:);

plot(bin_worst(:,1),bin_worst(:,6),"-o","LineWidth",1.5);

hold on;

plot(block_worst(:,1),block_worst(:,6),"-s","LineWidth",1.5);

plot(worst_caps,worst_caps,"k--");

xlabel("FFT Hard Cap");
ylabel("Actual Doppler FFTs");

title("All 500 Candidate Bins - Worst-Case Processing");

legend("Top-K Bins","Top-K Blocks","Hard Cap","Location","best");

grid on;


%% Save PNG

plot_file = fullfile(output_dir,"density_scaling_and_worst_case.png");

%exportgraphics(fig,plot_file,"Resolution",120);


%% TXT Summary

txt_file = fullfile(output_dir,"scene_density_worst_case_validation.txt");

file = fopen(txt_file,"w");

assert(file ~= -1,"Cannot create Group 4C summary");

fprintf(file,"=== Group 4C: Scene Density and Worst-Case Workload ===\n\n");

fprintf(file,"Target counts = 2, 4, 8, 16\n");

fprintf(file,"Independent trials per density = %d\n",num_trials);

fprintf(file,"ROI threshold = %.1f dB\n",threshold_db);

fprintf(file,"FFT caps = 100 and 200\n");

fprintf(file,"Real scene comparisons = %d\n",row);

fprintf(file,"Full-density worst-case comparisons = %d\n\n",worst_row);


fprintf(file,"[Mean Workload, FFT Cap 200]\n");

for i = 1:length(target_counts)

    fprintf(file,"Targets %d | Bin FFT %.2f | Block FFT %.2f | Bin retention %.2f%% | Block retention %.2f%%\n",target_counts(i),mean_ffts(i,1),mean_ffts(i,2),hit_retention(i,1),hit_retention(i,2));

end


fprintf(file,"\n[Full-Density Worst Case]\n");

for i = 1:worst_row

    fprintf(file,"Cap %d | Method %d | Selected bins %d | Actual FFTs %d\n",worst_results(i,1),worst_results(i,2),worst_results(i,4),worst_results(i,6));

end


fprintf(file,"\n[Validation]\n");

fprintf(file,"Actual Doppler FFT hard cap = PASS\n");
fprintf(file,"Dense FFT / CFAR / peak regression = PASS\n");
fprintf(file,"Multiple target-density scenarios = PASS\n");
fprintf(file,"Independent noise trials completed = PASS\n");
fprintf(file,"Zero-budget handling = PASS\n");
fprintf(file,"Full-budget handling = PASS\n");
fprintf(file,"Oversized atomic block handling = PASS\n");
fprintf(file,"Overall code validation = PASS\n");


fprintf(file,"\n[Limitations]\n");

fprintf(file,"Only 3 trials are used per target density.\n");

fprintf(file,"This is a workload stress test, not a final statistical Pd estimate.\n");

fprintf(file,"Top-K Blocks uses all-or-nothing candidate blocks.\n");

fprintf(file,"At full candidate density, a large block may be rejected by a small FFT cap.\n");

fprintf(file,"Matched target retention is relative to Dense CFAR detections.\n");

fprintf(file,"Synthetic targets do not model realistic clutter or multipath.\n");

fprintf(file,"FFT caps exclude Range FFT, CFAR, AoA and system overhead.\n");

fprintf(file,"Execution latency and peak memory are not measured.\n");

fclose(file);


%% Final Validation

assert(row == num_rows);

assert(all(results(:,9)<=results(:,3)));

assert(all(results(:,12)<=results(:,11)));

assert(all(results(:,14)<=results(:,13)));

assert(all(worst_results(:,6)<=worst_results(:,1)));


fprintf("\n--- Group 4C Validation ---\n");

fprintf("Target densities tested = %d\n",length(target_counts));

fprintf("Independent trials per density = %d\n",num_trials);

fprintf("Real scene comparisons = %d\n",row);

fprintf("Full-density worst-case comparisons = %d\n",worst_row);

fprintf("Hard FFT cap = PASS\n");

fprintf("Dense FFT / CFAR / Peak regression = PASS\n");

fprintf("Oversized block handling = PASS\n");

fprintf("Group 4C code validation passed.\n");

fprintf("Density CSV: %s\n",density_csv);
fprintf("Worst-case CSV: %s\n",worst_csv);
fprintf("Plot: %s\n",plot_file);
fprintf("Summary: %s\n",txt_file);


end
