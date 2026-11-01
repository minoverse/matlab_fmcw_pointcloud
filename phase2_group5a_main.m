
%% Phase 2 Group 5A - Unified Resource Counters

clear;
clc;

assert(isfile("config.m"), "Run from the project root");

run("config.m");


%% Input Result Files

file_4b = "output/phase2/09_topk_bounded_workload/topk_bounded_detection_results.csv";

file_4b2 = "output/phase2/10_topk_blocks_comparison/topk_bins_vs_blocks_comparison.csv";

file_4c = "output/phase2/11_scene_density_worst_case/scene_density_detection_workload.csv";

file_worst = "output/phase2/11_scene_density_worst_case/full_density_worst_case.csv";

assert(isfile(file_4b), "Group 4B CSV missing");
assert(isfile(file_4b2), "Group 4B-2 CSV missing");
assert(isfile(file_4c), "Group 4C density CSV missing");
assert(isfile(file_worst), "Group 4C worst-case CSV missing");


%% Read Previous Results

data_4b = readtable(file_4b);

data_4b2 = readtable(file_4b2);

data_4c = readtable(file_4c);

data_worst = readtable(file_worst);


%% Radar and CFAR Parameters

half_N = floor(Ns/2);

Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;

range_margin = Tr + Gr;
doppler_margin = Td + Gd;

training_cells = (2*range_margin+1)*(2*doppler_margin+1) - (2*Gr+1)*(2*Gd+1);

dense_cfar_cuts = (half_N-2*range_margin)*(Nchirp-2*doppler_margin);

current_dense_ffts = Ns*Nrx;
fair_dense_ffts = half_N*Nrx;

range_fft_count = Nrx*Nchirp;

assert(training_cells == 216);
assert(dense_cfar_cuts == 25920);
assert(current_dense_ffts == 4000);
assert(fair_dense_ffts == 2000);


%% FFT Workload Proxy

range_work = range_fft_count*Ns*log2(Ns);

doppler_work_per_fft = Nchirp*log2(Nchirp);

fair_dense_doppler_work = fair_dense_ffts*doppler_work_per_fft;

fair_total_fft_work = range_work + fair_dense_doppler_work;


%% Logical Data Sizes

complex_bytes = 16;
real_bytes = 8;
logical_bytes = 1;

adc_bytes = Ns*Nrx*Nchirp*complex_bytes;

range_cube_bytes = Ns*Nrx*Nchirp*complex_bytes;

range_fft_io_bytes = 2*range_fft_count*Ns*complex_bytes;

roi_energy_samples = half_N*Nrx*Nchirp;

roi_energy_read_bytes = roi_energy_samples*complex_bytes;


%% Build Unified Result Arrays

num_reference = 2;

num_rows = num_reference + height(data_4b) + height(data_4b2) + height(data_4c) + height(data_worst);

values = NaN(num_rows,12);

stage = strings(num_rows,1);
method = strings(num_rows,1);
profile = strings(num_rows,1);
hit_metric = strings(num_rows,1);

row = 0;


%% Reference 1: Original Dense Processing

row = row + 1;

stage(row) = "Reference";
method(row) = "Current Dense";
profile(row) = "Original";
hit_metric(row) = "Not recorded";

values(row,1) = current_dense_ffts;
values(row,2) = dense_cfar_cuts;
values(row,5) = Ns;


%% Reference 2: Fair Positive-Range Dense Processing

row = row + 1;

stage(row) = "Reference";
method(row) = "Fair Dense";
profile(row) = "Positive Range";
hit_metric(row) = "Not recorded";

values(row,1) = fair_dense_ffts;
values(row,2) = dense_cfar_cuts;
values(row,5) = half_N;


%% Group 4B Results

for i = 1:height(data_4b)

    row = row + 1;

    stage(row) = "Group 4B";
    method(row) = "Top-K Bins";
    profile(row) = string(data_4b.Profile(i));
    hit_metric(row) = "Both targets";

    values(row,:) = [
        data_4b.ActualDopplerFFTs(i), ...
        data_4b.CFARCUTs(i), ...
        data_4b.FFTCap(i), ...
        data_4b.TopKLimit(i), ...
        data_4b.SelectedROIBins(i), ...
        2, ...
        data_4b.Target2Amplitude(i), ...
        data_4b.ROIThreshold_dB(i), ...
        data_4b.DenseTarget1Hit(i)+data_4b.DenseTarget2Hit(i), ...
        data_4b.BoundedTarget1Hit(i)+data_4b.BoundedTarget2Hit(i), ...
        NaN, ...
        data_4b.CandidateROIBins(i)
    ];

end


%% Group 4B-2 Results

for i = 1:height(data_4b2)

    row = row + 1;

    stage(row) = "Group 4B-2";
    method(row) = string(data_4b2.Method(i));
    profile(row) = string(data_4b2.Profile(i));
    hit_metric(row) = "Target 2 only";

    values(row,:) = [
        data_4b2.ActualDopplerFFTs(i), ...
        data_4b2.CFARCUTs(i), ...
        data_4b2.FFTCap(i), ...
        data_4b2.TopKLimit(i), ...
        data_4b2.SelectedBins(i), ...
        2, ...
        data_4b2.Target2Amplitude(i), ...
        data_4b2.ROIThreshold_dB(i), ...
        data_4b2.DenseTarget2Hit(i), ...
        data_4b2.BoundedTarget2Hit(i), ...
        NaN, ...
        data_4b2.CandidateBins(i)
    ];

end


%% Group 4C Density Results

for i = 1:height(data_4c)

    row = row + 1;

    stage(row) = "Group 4C Density";
    method(row) = string(data_4c.Method(i));
    profile(row) = "Density Stress";
    hit_metric(row) = "All targets";

    values(row,:) = [
        data_4c.ActualDopplerFFTs(i), ...
        data_4c.CFARCUTs(i), ...
        data_4c.FFTCap(i), ...
        NaN, ...
        data_4c.SelectedBins(i), ...
        data_4c.TargetCount(i), ...
        NaN, ...
        NaN, ...
        data_4c.DenseTargetHits(i), ...
        data_4c.BoundedTargetHits(i), ...
        data_4c.Trial(i), ...
        data_4c.CandidateBins(i)
    ];

end


%% Group 4C Full-Density Worst Case

for i = 1:height(data_worst)

    row = row + 1;

    stage(row) = "Group 4C Worst";
    profile(row) = "Full Density";
    hit_metric(row) = "Not recorded";

    if data_worst.MethodIndex(i) == 1
        method(row) = "Bins";
    else
        method(row) = "Blocks";
    end

    values(row,:) = [
        data_worst.ActualDopplerFFTs(i), ...
        data_worst.CFARCUTs(i), ...
        data_worst.FFTCap(i), ...
        data_worst.TopKLimit(i), ...
        data_worst.SelectedBins(i), ...
        NaN, ...
        NaN, ...
        NaN, ...
        NaN, ...
        NaN, ...
        NaN, ...
        half_N
    ];

end


%% Validate Input Rows

assert(row == num_rows);

fft_counts = values(:,1);

cfar_cuts = values(:,2);

fft_caps = values(:,3);

selected_bins = values(:,5);

assert(all(isfinite(fft_counts)));
assert(all(isfinite(cfar_cuts)));

assert(all(fft_counts >= 0));
assert(all(cfar_cuts >= 0));

assert(all(mod(fft_counts,Nrx) == 0));

assert(all(cfar_cuts <= dense_cfar_cuts));

assert(all(fft_counts <= current_dense_ffts));

bounded_rows = 3:num_rows;

assert(all(fft_counts(bounded_rows) <= fft_caps(bounded_rows)));

assert(all(selected_bins(bounded_rows) <= half_N));

assert(all(values(bounded_rows,10) <= values(bounded_rows,9) | isnan(values(bounded_rows,9))));


%% Unified FFT Counters

expanded_bins = fft_counts/Nrx;

range_fft_counts = range_fft_count*ones(num_rows,1);

range_fft_work = range_work*ones(num_rows,1);

doppler_fft_work = fft_counts*doppler_work_per_fft;

total_fft_work = range_fft_work + doppler_fft_work;

doppler_reduction = 100*(1-fft_counts/fair_dense_ffts);

combined_fft_reduction = 100*(1-total_fft_work/fair_total_fft_work);


%% Unified CFAR Counters

cfar_training_accesses = cfar_cuts*training_cells;

cfar_reduction = 100*(1-cfar_cuts/dense_cfar_cuts);

cfar_logical_read_bytes = cfar_training_accesses*real_bytes;


%% Logical Memory Buffer Proxies

adc_mib = (adc_bytes/2^20)*ones(num_rows,1);

range_mib = (range_cube_bytes/2^20)*ones(num_rows,1);

doppler_output_bytes = fft_counts*Nchirp*complex_bytes;

doppler_output_mib = doppler_output_bytes/2^20;


%% RD Power Buffer

rd_processed_bins = expanded_bins;

rd_processed_bins(1) = half_N;
rd_processed_bins(2) = half_N;

rd_power_bytes = rd_processed_bins*Nchirp*real_bytes;

rd_power_mib = rd_power_bytes/2^20;


%% ROI CFAR Local Buffers

roi_noise_bytes = selected_bins*Nchirp*real_bytes;

roi_map_bytes = selected_bins*Nchirp*logical_bytes;

roi_noise_bytes(1:2) = half_N*Nchirp*real_bytes;

roi_map_bytes(1:2) = half_N*Nchirp*logical_bytes;

roi_local_mib = (roi_noise_bytes+roi_map_bytes)/2^20;


%% Logical Buffer Sum

buffer_sum_mib = adc_mib + range_mib + doppler_output_mib + rd_power_mib + roi_local_mib;


%% Logical FFT I/O

range_fft_io_mib = (range_fft_io_bytes/2^20)*ones(num_rows,1);

doppler_fft_io_bytes = 2*fft_counts*Nchirp*complex_bytes;

doppler_fft_io_mib = doppler_fft_io_bytes/2^20;


%% ROI Energy Scan

roi_scan_bytes = roi_energy_read_bytes*ones(num_rows,1);

roi_scan_bytes(1:2) = 0;

roi_scan_mib = roi_scan_bytes/2^20;


%% CFAR Logical Read Traffic

cfar_read_mib = cfar_logical_read_bytes/2^20;


%% Combined Logical Data Access Proxy

logical_data_io_mib = range_fft_io_mib + doppler_fft_io_mib + roi_scan_mib + cfar_read_mib;


%% Build Output Table

column_names = {
    'DopplerFFTs', ...
    'CFARCUTs', ...
    'FFTCap', ...
    'TopKLimit', ...
    'SelectedROIBins', ...
    'TargetCount', ...
    'Target2Amplitude', ...
    'ROIThreshold_dB', ...
    'ReportedDenseHits', ...
    'ReportedSparseHits', ...
    'Trial', ...
    'CandidateROIBins'
};

results_table = array2table(values,"VariableNames",column_names);

results_table.Stage = stage;
results_table.Method = method;
results_table.Profile = profile;
results_table.HitMetric = hit_metric;

results_table.RangeFFTCount = range_fft_counts;
results_table.RangeFFTWork = range_fft_work;

results_table.ExpandedRangeBins = expanded_bins;

results_table.DopplerFFTWork = doppler_fft_work;

results_table.TotalFFTWork = total_fft_work;

results_table.DopplerWorkReduction_pct = doppler_reduction;

results_table.TotalFFTWorkReduction_pct = combined_fft_reduction;

results_table.CFARTrainingAccesses = cfar_training_accesses;

results_table.CFARCUTReduction_pct = cfar_reduction;

results_table.ADCBuffer_MiB = adc_mib;
results_table.RangeBuffer_MiB = range_mib;
results_table.DopplerOutput_MiB = doppler_output_mib;
results_table.RDPowerBuffer_MiB = rd_power_mib;
results_table.CFARLocalBuffers_MiB = roi_local_mib;

results_table.LogicalBufferSum_MiB = buffer_sum_mib;

results_table.RangeFFTIO_MiB = range_fft_io_mib;
results_table.DopplerFFTIO_MiB = doppler_fft_io_mib;
results_table.ROIEnergyRead_MiB = roi_scan_mib;
results_table.CFARLogicalReads_MiB = cfar_read_mib;

results_table.LogicalDataAccess_MiB = logical_data_io_mib;


%% Output Folder

output_dir = "output/phase2/12_unified_resource_counters";

if ~isfolder(output_dir)
    mkdir(output_dir);
end


%% Save All Resource Counters

csv_file = fullfile(output_dir,"unified_resource_counters.csv");

writetable(results_table,csv_file);


%% Select Key Comparison Cases

reference_rows = stage == "Reference";

balanced_weak_rows = stage == "Group 4B-2" & profile == "Balanced" & abs(values(:,7)-0.005)<1e-12 & values(:,8)==3;

selected_rows = reference_rows | balanced_weak_rows;

comparison_table = results_table(selected_rows,:);

comparison_csv = fullfile(output_dir,"key_resource_comparison.csv");

writetable(comparison_table,comparison_csv);


%% Key Results

fprintf("\n--- Group 5A Unified Resource Counters ---\n");

fprintf("Input Group 4B rows = %d\n",height(data_4b));
fprintf("Input Group 4B-2 rows = %d\n",height(data_4b2));
fprintf("Input Group 4C density rows = %d\n",height(data_4c));
fprintf("Input Group 4C worst-case rows = %d\n",height(data_worst));

fprintf("Reference rows = %d\n",num_reference);
fprintf("Total resource rows = %d\n",num_rows);

fprintf("\nFair Dense:\n");

fprintf("Range FFT work = %.0f\n",range_work);
fprintf("Doppler FFT work = %.0f\n",fair_dense_doppler_work);
fprintf("Total FFT work proxy = %.0f\n",fair_total_fft_work);
fprintf("CFAR training accesses = %.0f\n",dense_cfar_cuts*training_cells);


%% Balanced Weak Scene

weak_indices = find(balanced_weak_rows);

assert(length(weak_indices) == 2);

for i = 1:length(weak_indices)

    k = weak_indices(i);

    fprintf("\n%s Balanced, Weak Target 0.005, ROI 3 dB:\n",method(k));

    fprintf("Doppler FFT count = %.0f\n",fft_counts(k));

    fprintf("Doppler FFT reduction = %.2f%%\n",doppler_reduction(k));

    fprintf("Total FFT work reduction = %.2f%%\n",combined_fft_reduction(k));

    fprintf("CFAR CUT reduction = %.2f%%\n",cfar_reduction(k));

    fprintf("Doppler output buffer = %.4f MiB\n",doppler_output_mib(k));

    fprintf("Logical data access proxy = %.4f MiB\n",logical_data_io_mib(k));

end


%% Validation

assert(all(abs(doppler_fft_work-fft_counts*doppler_work_per_fft)<1e-8));

assert(all(cfar_training_accesses==cfar_cuts*training_cells));

assert(all(buffer_sum_mib>=adc_mib+range_mib));

assert(all(logical_data_io_mib>=range_fft_io_mib));

assert(all(doppler_reduction(bounded_rows)>=0));

assert(all(combined_fft_reduction(bounded_rows)>=0));


%% TXT Summary

txt_file = fullfile(output_dir,"unified_resource_validation.txt");

file = fopen(txt_file,"w");

assert(file ~= -1,"Cannot create Group 5A summary");

fprintf(file,"=== Group 5A: Unified Resource Counters ===\n\n");

fprintf(file,"[Workload Settings]\n");

fprintf(file,"Range FFT size = %d\n",Ns);
fprintf(file,"Doppler FFT size = %d\n",Nchirp);
fprintf(file,"RX channels = %d\n",Nrx);
fprintf(file,"DSP Range FFT count = %d\n",range_fft_count);
fprintf(file,"Fair Dense Doppler FFT count = %d\n",fair_dense_ffts);
fprintf(file,"CFAR training cells = %d\n",training_cells);
fprintf(file,"Fair Dense CFAR CUTs = %d\n",dense_cfar_cuts);

fprintf(file,"\n[Input Results]\n");

fprintf(file,"Group 4B configurations = %d\n",height(data_4b));
fprintf(file,"Group 4B-2 configurations = %d\n",height(data_4b2));
fprintf(file,"Group 4C density configurations = %d\n",height(data_4c));
fprintf(file,"Group 4C worst-case configurations = %d\n",height(data_worst));
fprintf(file,"Total reference and experiment rows = %d\n",num_rows);

fprintf(file,"\n[FFT Work Proxies]\n");

fprintf(file,"Fair Dense Range FFT work = %.0f\n",range_work);
fprintf(file,"Fair Dense Doppler FFT work = %.0f\n",fair_dense_doppler_work);
fprintf(file,"Fair Dense combined FFT work = %.0f\n",fair_total_fft_work);

fprintf(file,"\n[Balanced Weak-Target Comparison]\n");

for i = 1:length(weak_indices)

    k = weak_indices(i);

    fprintf(file,"%s | Doppler FFTs %.0f | Doppler reduction %.2f%% | Combined FFT reduction %.2f%% | CFAR CUT reduction %.2f%%\n",method(k),fft_counts(k),doppler_reduction(k),combined_fft_reduction(k),cfar_reduction(k));

end

fprintf(file,"\n[Validation]\n");

fprintf(file,"Previous result files loaded = PASS\n");
fprintf(file,"Range FFT counters = PASS\n");
fprintf(file,"Doppler FFT workload counters = PASS\n");
fprintf(file,"CFAR training access counters = PASS\n");
fprintf(file,"FFT budget consistency = PASS\n");
fprintf(file,"Logical buffer size counters = PASS\n");
fprintf(file,"Logical data access counters = PASS\n");
fprintf(file,"Overall code validation = PASS\n");

fprintf(file,"\n[Important Limitations]\n");

fprintf(file,"FFT work uses N*log2(N), not measured CPU cycles.\n");
fprintf(file,"Combined FFT work does not include CFAR or AoA operations.\n");
fprintf(file,"Range FFT count excludes the original diagnostic duplicate FFT.\n");
fprintf(file,"CFAR accesses are logical training-cell reads, not measured DRAM traffic.\n");
fprintf(file,"Logical buffer sum is not measured peak RAM usage.\n");
fprintf(file,"Sparse processing buffers are modeled using selected range bins.\n");
fprintf(file,"Reported hit metrics differ across source experiments; do not aggregate them directly.\n");
fprintf(file,"AoA cost is excluded because per-scene grouped peak counts are not in the input CSV files.\n");
fprintf(file,"Original CFAR alpha was used in the Group 4 experiments.\n");
fprintf(file,"No actual CPU latency, energy or power reduction is claimed.\n");

fclose(file);


%% Final Output

fprintf("\n--- Group 5A Validation ---\n");

fprintf("Unified resource rows = %d\n",num_rows);

fprintf("Range FFT / Doppler FFT counters = PASS\n");
fprintf("CFAR training accesses = PASS\n");
fprintf("FFT budget consistency = PASS\n");
fprintf("Logical memory and data movement = PASS\n");

fprintf("Group 5A code validation passed.\n");

fprintf("Results CSV: %s\n",csv_file);
fprintf("Comparison CSV: %s\n",comparison_csv);
fprintf("Summary: %s\n",txt_file);

