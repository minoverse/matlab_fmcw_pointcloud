function phase2_group5b_main()

%% Phase 2 Group 5B - Performance vs Cost Comparison

% 함수의 지역 변수 범위가 자동으로 초기화됨
clc;

assert(isfile("src/config/radar_config.m"), "Run from the project root");

resource_file = "output/phase2/12_unified_resource_counters/unified_resource_counters.csv";
detection_file = "output/phase2/10_topk_blocks_comparison/topk_bins_vs_blocks_comparison.csv";

assert(isfile(resource_file), "Run Group 5A first");
assert(isfile(detection_file), "Run Group 4B-2 first");

resources = readtable(resource_file);
detections = readtable(detection_file);

resource_stage = string(resources.Stage);
resource_method = string(resources.Method);
resource_profile = string(resources.Profile);

detection_method = string(detections.Method);
detection_profile = string(detections.Profile);

reference_index = find(resource_stage == "Reference" & resource_method == "Fair Dense");

assert(length(reference_index) == 1);


%% Comparison Settings

amplitudes = [0.7, 0.005];
thresholds = [3, 6];

profiles = ["Strict", "Balanced", "Relaxed"];
methods = ["Bins", "Blocks"];

num_rows = length(amplitudes)*length(thresholds)*8;


%% Result Arrays

Method = strings(num_rows,1);
Profile = strings(num_rows,1);

Target2Amplitude = zeros(num_rows,1);
ROIThreshold_dB = zeros(num_rows,1);

FFTCap = zeros(num_rows,1);
DopplerFFTs = zeros(num_rows,1);
CFARCUTs = zeros(num_rows,1);
CFARTrainingAccesses = zeros(num_rows,1);

DopplerReduction_pct = zeros(num_rows,1);
TotalFFTReduction_pct = zeros(num_rows,1);

LogicalBufferSum_MiB = zeros(num_rows,1);
LogicalDataAccess_MiB = zeros(num_rows,1);

DenseTarget2Hit = zeros(num_rows,1);
SelectedTarget2Hit = zeros(num_rows,1);

Pareto = false(num_rows,1);

row = 0;


%% Paired Scene Comparison

for a = 1:length(amplitudes)

    for t = 1:length(thresholds)

        amplitude = amplitudes(a);
        threshold_db = thresholds(t);

        scenario = find(abs(detections.Target2Amplitude-amplitude)<1e-12 & detections.ROIThreshold_dB==threshold_db);

        assert(length(scenario)==8);

        sparse_source = find(abs(detections.Target2Amplitude-amplitude)<1e-12 & detections.ROIThreshold_dB==threshold_db & detection_profile=="Unbounded" & detection_method=="Bins");

        assert(length(sparse_source)==1);


        %% Dense, Sparse ROI, Top-K Bins, Top-K Blocks

        for choice = 1:8

            if choice == 1

                method_name = "Dense";
                profile_name = "Fair Dense";

                source_index = sparse_source;
                resource_index = reference_index;

                fft_cap = 2000;
                hit = detections.DenseTarget2Hit(source_index);

            elseif choice == 2

                method_name = "Sparse ROI";
                profile_name = "Unbounded";

                source_index = sparse_source;

                fft_cap = detections.FFTCap(source_index);
                hit = detections.BoundedTarget2Hit(source_index);

            else

                profile_index = floor((choice-3)/2)+1;
                method_index = mod(choice-3,2)+1;

                method_name = methods(method_index);
                profile_name = profiles(profile_index);

                source_index = find(abs(detections.Target2Amplitude-amplitude)<1e-12 & detections.ROIThreshold_dB==threshold_db & detection_profile==profile_name & detection_method==method_name);

                assert(length(source_index)==1);

                fft_cap = detections.FFTCap(source_index);
                hit = detections.BoundedTarget2Hit(source_index);

            end


            %% Match Group 5A Resource Data

            if choice ~= 1

                resource_index = find(resource_stage=="Group 4B-2" & resource_method==detection_method(source_index) & resource_profile==detection_profile(source_index) & abs(resources.Target2Amplitude-amplitude)<1e-12 & resources.ROIThreshold_dB==threshold_db);

                assert(length(resource_index)==1);

                assert(resources.DopplerFFTs(resource_index)==detections.ActualDopplerFFTs(source_index));

                assert(resources.CFARCUTs(resource_index)==detections.CFARCUTs(source_index));

            end


            %% Store Results

            row = row+1;

            Method(row) = method_name;
            Profile(row) = profile_name;

            Target2Amplitude(row) = amplitude;
            ROIThreshold_dB(row) = threshold_db;

            FFTCap(row) = fft_cap;

            DopplerFFTs(row) = resources.DopplerFFTs(resource_index);
            CFARCUTs(row) = resources.CFARCUTs(resource_index);

            CFARTrainingAccesses(row) = resources.CFARTrainingAccesses(resource_index);

            DopplerReduction_pct(row) = resources.DopplerWorkReduction_pct(resource_index);

            TotalFFTReduction_pct(row) = resources.TotalFFTWorkReduction_pct(resource_index);

            LogicalBufferSum_MiB(row) = resources.LogicalBufferSum_MiB(resource_index);

            LogicalDataAccess_MiB(row) = resources.LogicalDataAccess_MiB(resource_index);

            DenseTarget2Hit(row) = detections.DenseTarget2Hit(source_index);

            SelectedTarget2Hit(row) = hit;


            %% Validation

            assert(SelectedTarget2Hit(row)<=DenseTarget2Hit(row));

            assert(DopplerFFTs(row)<=FFTCap(row));

        end

    end

end

assert(row==num_rows);


%% Pareto Frontier Within Each Paired Scene

for i = 1:num_rows

    same_scene = abs(Target2Amplitude-Target2Amplitude(i))<1e-12 & ROIThreshold_dB==ROIThreshold_dB(i);

    dominated = false;

    candidates = find(same_scene);

    for j = 1:length(candidates)

        k = candidates(j);

        better_or_equal = DopplerFFTs(k)<=DopplerFFTs(i) & SelectedTarget2Hit(k)>=SelectedTarget2Hit(i);

        strictly_better = DopplerFFTs(k)<DopplerFFTs(i) | SelectedTarget2Hit(k)>SelectedTarget2Hit(i);

        if better_or_equal && strictly_better

            dominated = true;
            break;

        end

    end

    Pareto(i) = ~dominated;

end


%% Results Table

result_table = table(Method,Profile,Target2Amplitude,ROIThreshold_dB,FFTCap,DopplerFFTs,CFARCUTs,CFARTrainingAccesses,DopplerReduction_pct,TotalFFTReduction_pct,LogicalBufferSum_MiB,LogicalDataAccess_MiB,DenseTarget2Hit,SelectedTarget2Hit,Pareto);


%% Output Folder

output_dir = "output/phase2/13_performance_cost_tradeoff";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

results_csv = fullfile(output_dir,"performance_cost_comparison.csv");
pareto_csv = fullfile(output_dir,"pareto_frontier.csv");
summary_txt = fullfile(output_dir,"performance_cost_validation.txt");

writetable(result_table,results_csv);
writetable(result_table(Pareto,:),pareto_csv);


%% TXT Summary

file = fopen(summary_txt,"w");

assert(file~=-1,"Cannot save Group 5B summary");

fprintf(file,"=== Group 5B: Performance-Cost Tradeoff ===\n\n");

fprintf(file,"Paired scenario configurations = %d\n",num_rows);
fprintf(file,"Source = Group 4B-2 detection results and Group 5A resource counters\n");
fprintf(file,"CFAR uses original alpha, not calibrated alpha\n");

fprintf(file,"\n[Weak Target: Amplitude 0.005, ROI 3 dB]\n");

focus_rows = find(abs(Target2Amplitude-0.005)<1e-12 & ROIThreshold_dB==3);

for i = 1:length(focus_rows)

    k = focus_rows(i);

    fprintf(file,"%s / %s | FFT %d | Total FFT reduction %.2f%% | T2 hit %d/%d | Pareto %d\n",Method(k),Profile(k),DopplerFFTs(k),TotalFFTReduction_pct(k),SelectedTarget2Hit(k),DenseTarget2Hit(k),Pareto(k));

end

fprintf(file,"\n[Validation]\n");

fprintf(file,"Matched source rows = PASS\n");
fprintf(file,"Dense/Sparse paired-scene comparison = PASS\n");
fprintf(file,"Hard FFT budgets = PASS\n");
fprintf(file,"Pareto dominance comparison = PASS\n");
fprintf(file,"Overall code validation = PASS\n");

fprintf(file,"\n[Limitations]\n");

fprintf(file,"Each 4B-2 scene has one noise realization; hits are binary, not Pd estimates.\n");
fprintf(file,"Top-K limits count bins in Bins mode and source blocks in Blocks mode.\n");
fprintf(file,"Dense rows are paired references, not extra measured runs.\n");
fprintf(file,"FFT work is a mathematical proxy, not CPU cycles or execution time.\n");
fprintf(file,"Logical buffer sum and data access are not measured physical RAM or bus traffic.\n");
fprintf(file,"No PNG is generated to avoid MATLAB Online export failures.\n");

fclose(file);


%% Final Results

fprintf("\n--- Group 5B Validation ---\n");

fprintf("Paired comparison rows = %d\n",num_rows);

fprintf("Pareto rows = %d\n",nnz(Pareto));

fprintf("Source rows / budgets / paired comparisons = PASS\n");

fprintf("Group 5B code validation passed.\n");

fprintf("Results CSV: %s\n",results_csv);
fprintf("Pareto CSV: %s\n",pareto_csv);
fprintf("Summary: %s\n",summary_txt);

end
