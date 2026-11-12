function state = phase2_main()
clc;
close all;

%% 모듈 경로 등록

addpath(genpath("src"));


%% Radar 설정

cfg = radar_config();

c = cfg.c;
fc = cfg.fc;
B = cfg.B;
Tchirp = cfg.Tchirp;
fs = cfg.fs;
Nchirp = cfg.Nchirp;

Ns = cfg.Ns;
half_N = floor(Ns / 2);

S = cfg.S;
lambda = cfg.lambda;

Nrx = cfg.Nrx;
Nrx_az = cfg.Nrx_az;
Nrx_el = cfg.Nrx_el;

d_az = cfg.d_az;
d_el = cfg.d_el;

rx_y = cfg.rx_y;
rx_z = cfg.rx_z;

NangleFFT_az = cfg.NangleFFT_az;
NangleFFT_el = cfg.NangleFFT_el;

target_range = cfg.target_range;
target_velocity = cfg.target_velocity;
target_angle = cfg.target_angle;
target_elevation = cfg.target_elevation;
target_amplitude = cfg.target_amplitude;

Ntargets = cfg.Ntargets;

range_resolution = cfg.range_resolution;
tau = cfg.tau;
fb = cfg.fb;
fd = cfg.fd;

Tr = cfg.cfar.Tr;
Td = cfg.cfar.Td;
Gr = cfg.cfar.Gr;
Gd = cfg.cfar.Gd;
Pfa = cfg.cfar.Pfa;


%% ADC 신호 생성

snr_db = 20;

[adc_cube, t_fast] = simulate_adc_cube_process(cfg, snr_db);

plot_first_chirp_iq(adc_cube, t_fast);


%% Range FFT

[range_data, range_axis, range_data_no_window] = range_fft_process(adc_cube, cfg);

plot_range_fft(range_data, range_data_no_window, range_axis);


%% Pre-Doppler ROI 선택

roi_threshold_db = 6;

roi = roi_candidate_selection_process(range_data, half_N, range_axis, target_range, roi_threshold_db);

range_positive = range_data(1:half_N, :, :);

range_energy = roi.range_energy;
roi_noise_floor = roi.roi_noise_floor;
roi_threshold = roi.roi_threshold;

roi_mask = roi.roi_mask;
roi_bins = roi.roi_bins;

roi_count = roi.roi_count;
roi_fraction = roi.roi_fraction;

target_roi_selected = roi.target_roi_selected;

print_roi_candidate_selection(roi, half_N);

plot_roi_candidate_selection(roi, range_axis, "output/phase2/01_baseline_roi_gating");


%% 정적 클러터 처리

apply_clutter_removal = false;

range_data_for_doppler = static_clutter_process(range_data);


%% Dense Doppler FFT

[range_doppler_rx, RD_power, velocity_axis] = doppler_fft_process(range_data_for_doppler, cfg, apply_clutter_removal);

% 기존 Group 실험에서 사용하는 RD 변수 복원

RD_mag = sqrt(RD_power);

RD_db = 20*log10(RD_mag / max(RD_mag(:)) + eps);

plot_doppler_fft(RD_power, range_axis, velocity_axis, cfg);


%% Dense CA-CFAR

[cfar_map, noise_map, det_range_idx, det_doppler_idx, cfar_stats] = cfar2d_process(RD_power, cfg.cfar);

Ntraining = cfar_stats.Ntraining;
alpha = cfar_stats.alpha;

total_cells = cfar_stats.total_cells;
guard_cells = cfar_stats.guard_cells;

fprintf("\n--- 2D CA-CFAR ---\n");
fprintf("Training cells = %d\n", Ntraining);
fprintf("Pfa = %.1e\n", Pfa);
fprintf("Detected cells = %d\n", length(det_range_idx));

plot_cfar_detections(RD_power, det_range_idx, det_doppler_idx, range_axis, velocity_axis);


%% Peak Grouping

[peak_range_idx, peak_doppler_idx, peak_power] = peak_grouping_process(RD_power, det_range_idx, det_doppler_idx);

Npeaks = length(peak_range_idx);

fprintf("\n--- Peak Grouping ---\n");
fprintf("CFAR cells = %d\n", length(det_range_idx));
fprintf("Grouped peaks = %d\n", Npeaks);

for i = 1:Npeaks

    fprintf("Peak %d: Range = %.2f m, Velocity = %.2f m/s\n", i, range_axis(peak_range_idx(i)), velocity_axis(peak_doppler_idx(i)));

end


%% AoA 추정

[estimated_azimuth, estimated_elevation, estimated_snr] = aoa_fft_process(range_doppler_rx, RD_power, noise_map, peak_range_idx, peak_doppler_idx, cfg);

fprintf("\n--- Azimuth / Elevation AoA ---\n");

for i = 1:Npeaks

    fprintf("Detection %d: Azimuth = %.2f deg, Elevation = %.2f deg, SNR = %.2f dB\n", i, estimated_azimuth(i), estimated_elevation(i), estimated_snr(i));

end

plot_aoa_detections(estimated_azimuth, estimated_elevation);


%% Baseline Cost 계산

data.adc_cube = adc_cube;
data.range_data = range_data;
data.range_doppler_rx = range_doppler_rx;
data.RD_power = RD_power;

data.roi_count = roi_count;
data.roi_bins = roi_bins;

data.Npeaks = Npeaks;
data.peak_range_idx = peak_range_idx;

cost = baseline_cost_process(cfg, data);

print_baseline_cost(cost);


%% 기존 실험과의 호환 변수

range_fft_count = cost.range_fft_count;
range_fft_diagnostic_count = cost.range_fft_diagnostic_count;

doppler_fft_current = cost.doppler_fft_current;
doppler_fft_fair = cost.doppler_fft_fair;
doppler_fft_sparse = cost.doppler_fft_sparse;

projected_roi_reduction = cost.projected_roi_reduction;

cfar_cut_count = cost.cfar_cut_count;
cfar_training_accesses = cost.cfar_training_accesses;


%% 최종 검증

assert(all(target_roi_selected));
assert(roi_count < half_N);

assert(range_fft_count == Nrx * Nchirp);

assert(doppler_fft_current == Ns * Nrx);
assert(doppler_fft_fair == half_N * Nrx);
assert(doppler_fft_sparse == roi_count * Nrx);

assert(Ntraining == cost.Ntraining);
assert(cfar_training_accesses == cfar_cut_count * Ntraining);


%% Phase 2 검증 결과

fprintf("\n--- Phase 2 Validation ---\n");

fprintf("Target ROI validation = PASS\n");
fprintf("ROI reduction validation = PASS\n");
fprintf("Baseline counter validation = PASS\n");

fprintf("Selected range bins = %d / %d\n", roi_count, half_N);

fprintf("Projected ROI Doppler work reduction = %.2f%%\n", projected_roi_reduction);

fprintf("Phase 2 initial validation passed.\n");


%% 검증 결과 저장

output_dir = "output/phase2/01_baseline_roi_gating";

if ~isfolder(output_dir)
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
%% Export values for independent experiment functions
state = struct();
state.cfg = cfg;
state.c = c;
state.fc = fc;
state.B = B;
state.Tchirp = Tchirp;
state.fs = fs;
state.Nchirp = Nchirp;
state.Ns = Ns;
state.half_N = half_N;
state.S = S;
state.lambda = lambda;
state.Nrx = Nrx;
state.Nrx_az = Nrx_az;
state.Nrx_el = Nrx_el;
state.d_az = d_az;
state.d_el = d_el;
state.rx_y = rx_y;
state.rx_z = rx_z;
state.NangleFFT_az = NangleFFT_az;
state.NangleFFT_el = NangleFFT_el;
state.target_range = target_range;
state.target_velocity = target_velocity;
state.target_angle = target_angle;
state.target_elevation = target_elevation;
state.target_amplitude = target_amplitude;
state.Ntargets = Ntargets;
state.range_resolution = range_resolution;
state.tau = tau;
state.fb = fb;
state.fd = fd;
state.Tr = Tr;
state.Td = Td;
state.Gr = Gr;
state.Gd = Gd;
state.Pfa = Pfa;
state.snr_db = snr_db;
state.adc_cube = adc_cube;
state.t_fast = t_fast;
state.range_data = range_data;
state.range_axis = range_axis;
state.range_data_no_window = range_data_no_window;
state.roi_threshold_db = roi_threshold_db;
state.roi = roi;
state.range_positive = range_positive;
state.range_energy = range_energy;
state.roi_noise_floor = roi_noise_floor;
state.roi_threshold = roi_threshold;
state.roi_mask = roi_mask;
state.roi_bins = roi_bins;
state.roi_count = roi_count;
state.roi_fraction = roi_fraction;
state.target_roi_selected = target_roi_selected;
state.apply_clutter_removal = apply_clutter_removal;
state.range_data_for_doppler = range_data_for_doppler;
state.range_doppler_rx = range_doppler_rx;
state.RD_power = RD_power;
state.velocity_axis = velocity_axis;
state.RD_mag = RD_mag;
state.RD_db = RD_db;
state.cfar_map = cfar_map;
state.noise_map = noise_map;
state.det_range_idx = det_range_idx;
state.det_doppler_idx = det_doppler_idx;
state.cfar_stats = cfar_stats;
state.Ntraining = Ntraining;
state.alpha = alpha;
state.total_cells = total_cells;
state.guard_cells = guard_cells;
state.peak_range_idx = peak_range_idx;
state.peak_doppler_idx = peak_doppler_idx;
state.peak_power = peak_power;
state.Npeaks = Npeaks;
state.estimated_azimuth = estimated_azimuth;
state.estimated_elevation = estimated_elevation;
state.estimated_snr = estimated_snr;
state.data = data;
state.cost = cost;
state.range_fft_count = range_fft_count;
state.range_fft_diagnostic_count = range_fft_diagnostic_count;
state.doppler_fft_current = doppler_fft_current;
state.doppler_fft_fair = doppler_fft_fair;
state.doppler_fft_sparse = doppler_fft_sparse;
state.projected_roi_reduction = projected_roi_reduction;
state.cfar_cut_count = cfar_cut_count;
state.cfar_training_accesses = cfar_training_accesses;
end
