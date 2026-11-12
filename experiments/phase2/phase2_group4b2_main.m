function phase2_group4b2_main()

%% Phase 2 Group 4B-2 - Top-K Blocks Comparison

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

half_N = radar.half_N;
Tr = radar.Tr;
Gr = radar.Gr;


%% Synthetic Unit Tests

roi_block_budget_regression(half_N,Nrx,Tr,Gr);


%% Experiment Settings

scene_amplitudes = [0.7,0.005];
roi_thresholds = [3,6];

profile_names = ["Strict","Balanced","Relaxed","Unbounded"];

top_k_limits = [3,6,12,half_N];
fft_caps = [100,200,400,half_N*Nrx];

methods = ["Bins","Blocks"];

num_rows = length(scene_amplitudes)*length(roi_thresholds)*length(profile_names)*length(methods);

results = zeros(num_rows,18);

method_labels = strings(num_rows,1);
profile_labels = strings(num_rows,1);

row = 0;

rng(12);


%% Real Scene Comparison

for scene = 1:length(scene_amplitudes)

    amplitudes = [1,scene_amplitudes(scene)];

    data = roi_budget_scene_data(target_range,target_velocity,amplitudes,radar);

    range_energy = data.range_energy;

    noise_floor = median(range_energy);

    fprintf("\n--- Group 4B-2 Scene: amplitude %.3f ---\n",amplitudes(2));


    for t = 1:length(roi_thresholds)

        threshold_db = roi_thresholds(t);

        threshold = noise_floor*10^(threshold_db/10);

        candidate_bins = find(range_energy > threshold);


        for p = 1:length(profile_names)

            top_k = top_k_limits(p);
            fft_cap = fft_caps(p);


            for m = 1:length(methods)

                method = methods(m);

                result = roi_budget_compare_case(method,candidate_bins,data,radar,top_k,fft_cap);

                row = row+1;

                method_labels(row) = method;
                profile_labels(row) = profile_names(p);

                results(row,:) = [amplitudes(2),threshold_db,p,top_k,fft_cap,result.candidate_bins,result.candidate_blocks,m,result.selected_units,result.selected_bins,result.selected_source_blocks,result.processing_blocks,result.fft_count,result.cfar_cuts,result.dense_hits(end),result.bounded_hits(end),result.fft_reduction,result.dropped_bins];

                fprintf("%s | %s | ROI %.1f dB | FFT %d/%d | Target 2 %d/%d\n",method,profile_names(p),threshold_db,result.fft_count,fft_cap,result.bounded_hits(end),result.dense_hits(end));

            end

        end

    end

end


%% Results Table

assert(row == num_rows);

column_names = {
    'Target2Amplitude','ROIThreshold_dB','ProfileIndex','TopKLimit','FFTCap','CandidateBins','CandidateBlocks','MethodIndex','SelectedUnits','SelectedBins','SelectedSourceBlocks','ProcessingBlocks','ActualDopplerFFTs','CFARCUTs','DenseTarget2Hit','BoundedTarget2Hit','FFTWorkloadReduction_pct','DroppedCandidateBins'
};

results_table = array2table(results,"VariableNames",column_names);

results_table.Method = method_labels;
results_table.Profile = profile_labels;


%% Output Folder

output_dir = "output/phase2/10_topk_blocks_comparison";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

csv_file = fullfile(output_dir,"topk_bins_vs_blocks_comparison.csv");

writetable(results_table,csv_file);


%% Weak Scene Comparison Plot

weak_bin_rows = find(results(:,1)==0.005 & results(:,2)==3 & results(:,8)==1);

weak_block_rows = find(results(:,1)==0.005 & results(:,2)==3 & results(:,8)==2);

assert(length(weak_bin_rows)==length(profile_names));

assert(length(weak_block_rows)==length(profile_names));


fig = figure("Position",[100 100 1000 700]);


%% Doppler FFT Workload

subplot(2,1,1);

plot(1:length(profile_names),results(weak_bin_rows,13),"-o","LineWidth",1.5);

hold on;

plot(1:length(profile_names),results(weak_block_rows,13),"-s","LineWidth",1.5);

plot(1:length(profile_names),fft_caps,"k--");

xticks(1:length(profile_names));
xticklabels(profile_names);

ylabel("Doppler FFT Count");

title("Top-K Bins vs Top-K Blocks: Same FFT Budget");

legend("Top-K Bins","Top-K Blocks","Hard FFT Cap","Location","best");

grid on;


%% Weak Target Detection

subplot(2,1,2);

plot(1:length(profile_names),results(weak_bin_rows,16),"-o","LineWidth",1.5);

hold on;

plot(1:length(profile_names),results(weak_block_rows,16),"-s","LineWidth",1.5);

plot(1:length(profile_names),results(weak_bin_rows,15),"k--");

xticks(1:length(profile_names));
xticklabels(profile_names);

ylabel("Weak Target Hit");

xlabel("Workload Profile");

title("Weak Target Preservation: Amplitude 0.005, ROI 3 dB");

legend("Top-K Bins","Top-K Blocks","Dense CFAR","Location","best");

ylim([-0.1 1.1]);

grid on;


%% Save Low-Resolution PNG

plot_file = fullfile(output_dir,"topk_bins_vs_blocks_tradeoff.png");

exportgraphics(fig,plot_file,"Resolution",120);


%% TXT Summary

txt_file = fullfile(output_dir,"topk_blocks_comparison_validation.txt");

file = fopen(txt_file,"w");

assert(file ~= -1,"Cannot create Group 4B-2 summary");

fprintf(file,"=== Group 4B-2: Top-K Bins vs Top-K Blocks ===\n\n");

fprintf(file,"Real scene configurations = %d\n",row);
fprintf(file,"CFAR range margin = %d\n",Tr+Gr);
fprintf(file,"Fair Dense FFTs = %d\n",half_N*Nrx);
fprintf(file,"CFAR alpha = %.6f\n",radar.Ntraining*(radar.Pfa^(-1/radar.Ntraining)-1));

fprintf(file,"\n[Weak Target: Amplitude 0.005, ROI 3 dB]\n");

for p = 1:length(profile_names)

    a = weak_bin_rows(p);
    b = weak_block_rows(p);

    fprintf(file,"%s | Cap %d | Bin FFT %d, Hit %d | Block FFT %d, Hit %d\n",profile_names(p),fft_caps(p),results(a,13),results(a,16),results(b,13),results(b,16));

end

fprintf(file,"\n[Validation]\n");

fprintf(file,"Top-K Bins regression = PASS\n");
fprintf(file,"Top-K Blocks regression = PASS\n");
fprintf(file,"Hard FFT budget = PASS\n");
fprintf(file,"CFAR margin coverage = PASS\n");
fprintf(file,"Dense FFT / CFAR / peak regression = PASS\n");
fprintf(file,"Synthetic weak-target coverage test = PASS\n");
fprintf(file,"Oversized block handling = PASS\n");
fprintf(file,"Overall code validation = PASS\n");

fprintf(file,"\n[Limitations]\n");

fprintf(file,"Top-K counts different units for Bins and Blocks.\n");
fprintf(file,"Both methods share the same FFT cap.\n");
fprintf(file,"Block score is maximum candidate range energy.\n");
fprintf(file,"Block selection is greedy and does not guarantee target preservation.\n");
fprintf(file,"Source blocks are not split when they exceed the FFT budget.\n");
fprintf(file,"Each real scene uses one independent noise realization.\n");
fprintf(file,"FFT reduction is not a measured latency or power saving.\n");

fclose(file);


%% Final Validation

assert(all(results(:,13)<=results(:,5)));

assert(all(results(:,9)<=results(:,4)));

assert(all(results(:,16)<=results(:,15)));

fprintf("\n--- Group 4B-2 Validation ---\n");

fprintf("Real scene comparisons = %d\n",row);

fprintf("Top-K Bins vs Blocks = PASS\n");
fprintf("Hard FFT cap = PASS\n");
fprintf("Dense FFT / CFAR / Peak regression = PASS\n");
fprintf("Group 4B-2 code validation passed.\n");

fprintf("CSV: %s\n",csv_file);
fprintf("PNG: %s\n",plot_file);
fprintf("TXT: %s\n",txt_file);

end
