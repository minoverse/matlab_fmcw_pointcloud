
clear;
clc;
close all;

addpath(genpath("src"));

run("config.m");


%% Simulation

run("simulate_adc_cube.m");


%% Range FFT

run("range_fft.m");


%% Pre-Doppler ROI Gating

run("roi_candidate_selection.m");


%% Dense Static / DC Handling

run("static_clutter_removal.m");


%% Dense Doppler FFT

run("doppler_fft.m");


%% Detection

run("cfar2d.m");

run("peak_grouping.m");


%% Angle Estimation

run("aoa_fft.m");


%% Baseline Cost Breakdown

run("baseline_cost_breakdown.m");


%% Final Validation

assert(all(target_roi_selected));
assert(roi_count < half_N);

assert(range_fft_count == Nrx * Nchirp);

assert(doppler_fft_current == Ns * Nrx);
assert(doppler_fft_fair == half_N * Nrx);
assert(doppler_fft_sparse == roi_count * Nrx);

assert(cfar_training_accesses == cfar_cut_count * Ntraining);


%% Results

fprintf("\n--- Phase 2 Validation ---\n");

fprintf("Target ROI validation = PASS\n");

fprintf("ROI reduction validation = PASS\n");

fprintf("Baseline counter validation = PASS\n");

fprintf("Selected range bins = %d / %d\n", roi_count, half_N);

fprintf("Projected ROI Doppler work reduction = %.2f%%\n", projected_roi_reduction);

fprintf("Phase 2 initial validation passed.\n");

%% Save Baseline Cost and ROI Validation Results

output_dir = "output/phase2/01_baseline_roi_gating";

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

result_file = fullfile(output_dir, "baseline_cost_roi_validation.txt");

file = fopen(result_file, "w");

assert(file ~= -1, "Cannot create result file");

fprintf(file, "=== Phase 2: Baseline Cost and Pre-Doppler ROI Validation ===\n\n");

fprintf(file, "[ROI Candidate Selection]\n");
fprintf(file, "Threshold = %.2f dB\n", roi_threshold_db);
fprintf(file, "Selected range bins = %d / %d\n", roi_count, half_N);
fprintf(file, "Selected fraction = %.2f%%\n", roi_fraction * 100);

for k = 1:Ntargets
    fprintf(file, "Target %d ROI selected = %d\n", k, target_roi_selected(k));
end

fprintf(file, "\n[Range FFT Baseline]\n");
fprintf(file, "DSP Range FFT count = %d\n", range_fft_count);
fprintf(file, "Diagnostic Range FFT count = %d\n", range_fft_diagnostic_count);

fprintf(file, "\n[Doppler FFT Workload]\n");
fprintf(file, "Current Dense FFT count = %d\n", doppler_fft_current);
fprintf(file, "Fair Dense FFT count = %d\n", doppler_fft_fair);
fprintf(file, "Projected Sparse FFT count = %d\n", doppler_fft_sparse);
fprintf(file, "Projected Doppler FFT work reduction = %.2f%%\n", projected_roi_reduction);

fprintf(file, "\n[CFAR Baseline]\n");
fprintf(file, "Training cells = %d\n", Ntraining);
fprintf(file, "CFAR CUT count = %d\n", cfar_cut_count);
fprintf(file, "Training accesses = %d\n", cfar_training_accesses);

fprintf(file, "\n[Validation]\n");
fprintf(file, "Target ROI validation = PASS\n");
fprintf(file, "ROI reduction validation = PASS\n");
fprintf(file, "Baseline counter validation = PASS\n");
fprintf(file, "Overall = PASS\n");

fprintf(file, "\n[Note]\n");
fprintf(file, "Sparse Doppler cost is a workload projection only.\n");
fprintf(file, "No measured CPU latency or power reduction is claimed.\n");

fclose(file);

fprintf("Validation results saved: %s\n", result_file);
