function phase2_group5d_main()
%% Phase 2 Group 5D - Fair CFAR Data Access Comparison
% New validation only. Existing signal processing functions are unchanged.
% Estimates logical accesses, not measured SRAM/DRAM traffic or latency.

% 함수의 지역 변수 범위가 자동으로 초기화됨
clc;

assert(isfile("src/config/radar_config.m"),"Run from the project root");
addpath(genpath("src"));
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

%% Radar / CFAR settings
half_N = floor(Ns/2);
Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;
margin_r = Tr+Gr;
margin_d = Td+Gd;
Ntraining = (2*margin_r+1)*(2*margin_d+1)-(2*Gr+1)*(2*Gd+1);
Pfa = 1e-5;
alpha = Ntraining*(Pfa^(-1/Ntraining)-1);
assert(Ntraining==216);

radar = struct();
radar.Ns = Ns;
radar.Nrx = Nrx;
radar.Nchirp = Nchirp;
radar.half_N = half_N;
radar.fs = fs;
radar.Tchirp = Tchirp;
radar.c = c;
radar.S = S;
radar.lambda = lambda;
radar.rx_y = rx_y;
radar.rx_z = rx_z;
radar.Tr = Tr;
radar.Td = Td;
radar.Gr = Gr;
radar.Gd = Gd;
radar.Ntraining = Ntraining;
radar.Pfa = Pfa;
radar.snr_db = 20;

%% Paired reference scene: target 2 = 0.005, ROI threshold 3 dB
rng(51005);
data = roi_budget_scene_data(target_range,target_velocity,[1,0.005],radar);
range_energy = data.range_energy;
threshold_db = 3;
roi_candidates = find(range_energy > median(range_energy)*10^(threshold_db/10));
[selected_bins,blocks,expanded_bins,~] = roi_budget_selection(roi_candidates,range_energy,half_N,margin_r,6,200,Nrx);
[sparse_fft,sparse_power,actual_fft_count] = block_doppler_processing(data.range_data,blocks,false);
assert(actual_fft_count==numel(expanded_bins)*Nrx);
assert(actual_fft_count<=200);

%% Validate actual sparse Doppler processing
if ~isempty(expanded_bins)
    expected_fft = data.dense_fft(expanded_bins,:,:);
    fft_error = max(abs(sparse_fft(:)-expected_fft(:)));
    assert(fft_error <= 1e-10*max(1,max(abs(expected_fft(:)))));
end

%% Dense CFAR with a 2-D summed-area table (integral image)
dense_power = data.RD_power;
P_dense = zeros(half_N+1,Nchirp+1);
P_dense(2:end,2:end) = cumsum(cumsum(dense_power,1),2);
[R,D] = ndgrid(margin_r+1:half_N-margin_r,margin_d+1:Nchirp-margin_d);
full_sum = rect_sum(P_dense,R-margin_r,R+margin_r,D-margin_d,D+margin_d);
guard_sum = rect_sum(P_dense,R-Gr,R+Gr,D-Gd,D+Gd);
dense_noise = (full_sum-guard_sum)/Ntraining;
valid_dense_noise = data.noise_map(R(:,1),D(1,:));
noise_error = max(abs(dense_noise(:)-valid_dense_noise(:)));
assert(noise_error <= 1e-8*max(1,max(abs(valid_dense_noise(:)))));
dense_detection = dense_power(R(:,1),D(1,:)) > alpha*dense_noise;
expected_dense_detection = data.cfar_map(R(:,1),D(1,:));
assert(isequal(dense_detection,expected_dense_detection));

%% Sparse integral-image CFAR, processed block by processed block
sparse_noise = zeros(numel(selected_bins),Nchirp);
sparse_detection = false(numel(selected_bins),Nchirp);
valid_d = margin_d+1:Nchirp-margin_d;
position = 1;

for b = 1:size(blocks,1)
    first_bin = blocks(b,1);
    last_bin = blocks(b,2);
    width = last_bin-first_bin+1;
    block_power = sparse_power(position:position+width-1,:);
    P_block = zeros(width+1,Nchirp+1);
    P_block(2:end,2:end) = cumsum(cumsum(block_power,1),2);
    idx = find(selected_bins>=first_bin & selected_bins<=last_bin & selected_bins>margin_r & selected_bins<=half_N-margin_r);

    for ii = 1:length(idx)
        original_r = selected_bins(idx(ii));
        local_r = original_r-first_bin+1;
        assert(local_r-margin_r>=1 && local_r+margin_r<=width);
        d = valid_d;
        full_values = rect_sum(P_block,local_r-margin_r,local_r+margin_r,d-margin_d,d+margin_d);
        guard_values = rect_sum(P_block,local_r-Gr,local_r+Gr,d-Gd,d+Gd);
        sparse_noise(idx(ii),d) = (full_values-guard_values)/Ntraining;
        sparse_detection(idx(ii),d) = block_power(local_r,d) > alpha*sparse_noise(idx(ii),d);
    end
    position = position+width;
end

assert(position==numel(expanded_bins)+1);
expected_sparse_noise = data.noise_map(selected_bins,:);
noise_difference = abs(sparse_noise-expected_sparse_noise);
if ~isempty(noise_difference)
    assert(max(noise_difference(:)) <= 1e-8*max(1,max(abs(expected_sparse_noise(:)))));
end
assert(isequal(sparse_detection,data.cfar_map(selected_bins,:)));

%% Fair logical I/O model: same CFAR algorithm on Dense and Sparse
% Prefix model: two cumulative-sum passes, each reads and writes every cell.
% CFAR query: 4 full-window corners + 4 guard-window corners + 1 CUT.
% Array allocation, cache behavior, burst effects and RAM bus transactions excluded.
complex_bytes = 16;
real_bytes = 8;
range_fft_io = 2*Nrx*Nchirp*Ns*complex_bytes;
dense_doppler_fft_io = 2*half_N*Nrx*Nchirp*complex_bytes;
sparse_doppler_fft_io = 2*actual_fft_count*Nchirp*complex_bytes;
roi_energy_read = half_N*Nrx*Nchirp*complex_bytes;
dense_cuts = (half_N-2*margin_r)*(Nchirp-2*margin_d);
valid_selected = selected_bins(selected_bins>margin_r & selected_bins<=half_N-margin_r);
sparse_cuts = numel(valid_selected)*(Nchirp-2*margin_d);

methods = ["Fair Dense";"Balanced Sparse"];
ffts = [half_N*Nrx;actual_fft_count];
processed_bins = [half_N;numel(expanded_bins)];
cuts = [dense_cuts;sparse_cuts];
range_io = [range_fft_io;range_fft_io];
doppler_io = [dense_doppler_fft_io;sparse_doppler_fft_io];
gating_io = [0;roi_energy_read];
naive_cfar_io = cuts*Ntraining*real_bytes;
prefix_build_io = 4*processed_bins*Nchirp*real_bytes;
prefix_query_io = 9*cuts*real_bytes;
prefix_cfar_io = prefix_build_io+prefix_query_io;
naive_total = range_io+doppler_io+gating_io+naive_cfar_io;
prefix_total = range_io+doppler_io+gating_io+prefix_cfar_io;

results = table(methods,ffts,processed_bins,cuts,range_io/2^20,doppler_io/2^20,gating_io/2^20,naive_cfar_io/2^20,prefix_build_io/2^20,prefix_query_io/2^20,prefix_cfar_io/2^20,naive_total/2^20,prefix_total/2^20);
results.Properties.VariableNames = {'Method','DopplerFFTs','ProcessedRangeBins','CFARCUTs','RangeFFTIO_MiB','DopplerFFTIO_MiB','ROIEnergyRead_MiB','NaiveCFARRead_MiB','PrefixCFARBuildIO_MiB','PrefixCFARQueryRead_MiB','PrefixCFARTotalIO_MiB','NaiveTotalAccess_MiB','PrefixTotalAccess_MiB'};

naive_saving = 100*(1-naive_total(2)/naive_total(1));
prefix_saving = 100*(1-prefix_total(2)/prefix_total(1));
assert(prefix_saving<naive_saving);
assert(all(prefix_cfar_io>=0));

%% Save results
output_dir = "output/phase2/15_fair_cfar_cost_validation";
if ~isfolder(output_dir)
    mkdir(output_dir);
end
csv_file = fullfile(output_dir,"fair_cfar_cost_comparison.csv");
txt_file = fullfile(output_dir,"fair_cfar_cost_validation.txt");
writetable(results,csv_file);

file = fopen(txt_file,"w");
assert(file~=-1,"Cannot write Group 5D TXT");
fprintf(file,"=== Group 5D: Fair CFAR Data Access Comparison ===\n\n");
fprintf(file,"Scene: SNR 20 dB, amplitude 0.005, ROI 3 dB, Top-K 6, cap 200\n");
fprintf(file,"Actual sparse Doppler FFTs = %d\n",actual_fft_count);
fprintf(file,"Selected CUTs = %d\n",sparse_cuts);
fprintf(file,"Naive CFAR model access reduction = %.2f%%\n",naive_saving);
fprintf(file,"Equal-method integral-image model access reduction = %.2f%%\n",prefix_saving);
fprintf(file,"Dense prefix-model logical access = %.4f MiB\n",prefix_total(1)/2^20);
fprintf(file,"Sparse prefix-model logical access = %.4f MiB\n",prefix_total(2)/2^20);
fprintf(file,"Dense/sparse FFT and integral CFAR regression = PASS\n\n");
fprintf(file,"[Limitations]\n");
fprintf(file,"Model: two prefix passes (read/write), 8 window-corner reads plus 1 CUT read.\n");
fprintf(file,"Both Dense and Sparse use the same integral-image cost model.\n");
fprintf(file,"Real CFAR kernels may use sliding sums, vectorization, caching, or tiling.\n");
fprintf(file,"ROI median, sorting, scheduling and block merge computation are not counted.\n");
fprintf(file,"MATLAB allocations and true physical bus bytes/latency are not measured.\n");
fprintf(file,"A modeled reduction is NOT an upper bound on all hardware implementations.\n");
fclose(file);

fprintf("\n--- Group 5D Validation ---\n");
fprintf("Naive CFAR model: %.2f%% reduction\n",naive_saving);
fprintf("Fair prefix CFAR model: %.2f%% reduction\n",prefix_saving);
fprintf("Actual sparse FFT / prefix CFAR regression = PASS\n");
fprintf("Group 5D code validation passed.\n");
fprintf("Results CSV: %s\n",csv_file);
fprintf("Summary: %s\n",txt_file);

%% Local rectangle-sum helper (1-based input rectangle coordinates)
end

function values = rect_sum(P,r1,r2,c1,c2)
    if isscalar(r1) && ~isscalar(c1)
        r1 = repmat(r1,size(c1));
        r2 = repmat(r2,size(c1));
    end
    values = P(sub2ind(size(P),r2+1,c2+1))-P(sub2ind(size(P),r1,c2+1))-P(sub2ind(size(P),r2+1,c1))+P(sub2ind(size(P),r1,c1));
end
