%% Phase 2 Group 5E - Dense vs Actual Sparse Range/Velocity/AoA RMSE
% Uses the same ADC scene for both paths; no figures or exportgraphics.
% The existing roi_only_cfar.m uses original alpha approximately 11.825, NOT 4.8.
% Reports both matched conditional RMSE and missed targets.

clear;
clc;

assert(isfile("config.m"),"Run from the project root");
addpath(genpath("src"));
run("config.m");
assert(Ntargets==2,"Group 5E uses the existing two-target geometry");

%% Experiment configuration
scene_names = ["20dB ROI3";"20dB ROI6";"10dB ROI3"];
scene_snr = [20,20,10];
scene_threshold = [3,6,3];
weak_amplitude = 0.005;
num_trials = 30;
fft_cap = 200;
top_k = 6;
num_scenes = length(scene_snr);
num_methods = 2;
num_metrics = 4;

%% Radar configuration
half_N = floor(Ns/2);
Tr = 8;
Td = 4;
Gr = 2;
Gd = 1;
margin = Tr+Gr;
Ntraining = (2*(Tr+Gr)+1)*(2*(Td+Gd)+1)-(2*Gr+1)*(2*Gd+1);
Pfa = 1e-5;
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

%% Target indices used for one-to-one matching within +/-1 bin
range_axis = c*(0:half_N-1)*fs/(2*S*Ns);
velocity_axis = (-Nchirp/2:Nchirp/2-1)*(1/Tchirp)*lambda/(2*Nchirp);
true_range_bins = zeros(Ntargets,1);
true_doppler_bins = zeros(Ntargets,1);
for k = 1:Ntargets
    [~,true_range_bins(k)] = min(abs(range_axis-target_range(k)));
    [~,true_doppler_bins(k)] = min(abs(velocity_axis-target_velocity(k)));
end

%% Statistics
hits = zeros(num_scenes,num_methods);
common_hits = zeros(num_scenes,1);
errors_all = zeros(num_scenes,num_methods,num_metrics);
errors_common = zeros(num_scenes,num_methods,num_metrics);
fft_total = zeros(num_scenes,num_methods);
cfar_total = zeros(num_scenes,num_methods);
max_estimate_difference = 0;

%% Paired Dense / Sparse actual processing
for s = 1:num_scenes
    radar.snr_db = scene_snr(s);
    % 20dB ROI3 and 20dB ROI6 use the same noise draws (paired comparison).
    rng(5500+scene_snr(s));

    for trial = 1:num_trials
        data = roi_budget_scene_data(target_range,target_velocity,[1,weak_amplitude],radar);

        %% Actual Sparse FFT with balanced Top-K cap
        candidate_bins = find(data.range_energy > median(data.range_energy)*10^(scene_threshold(s)/10));
        [selected_bins,processing_blocks,expanded_bins,~] = roi_budget_selection(candidate_bins,data.range_energy,half_N,margin,top_k,fft_cap,Nrx);
        [sparse_fft,sparse_RD_power,actual_ffts] = block_doppler_processing(data.range_data,processing_blocks,false);
        assert(actual_ffts <= fft_cap);
        assert(actual_ffts == numel(expanded_bins)*Nrx);

        if ~isempty(expanded_bins)
            expected_fft = data.dense_fft(expanded_bins,:,:);
            fft_difference = max(abs(sparse_fft(:)-expected_fft(:)));
            assert(fft_difference <= 1e-10*max(1,max(abs(expected_fft(:)))));
        end

        %% Actual ROI-only CFAR; original alpha from the repository
        roi_bins = selected_bins;
        active_range_bins = expanded_bins;
        evalc('run("roi_only_cfar.m");');
        assert(isequal(roi_cfar_local,data.cfar_map(selected_bins,:)));

        %% Dense local peak grouping (same 3x3 criterion)
        [dense_r,dense_d] = find(data.cfar_map);
        dense_peaks = zeros(0,2);
        for i = 1:length(dense_r)
            r = dense_r(i);
            d = dense_d(i);
            if r<=1 || r>=half_N || d<=1 || d>=Nchirp
                continue;
            end
            local_area = data.RD_power(r-1:r+1,d-1:d+1);
            if data.RD_power(r,d) == max(local_area(:))
                dense_peaks(end+1,:) = [r,d];
            end
        end
        sparse_peaks = [roi_peak_range_idx(:),roi_peak_doppler_idx(:)];

        %% Workload counters
        fft_total(s,1) = fft_total(s,1)+half_N*Nrx;
        fft_total(s,2) = fft_total(s,2)+actual_ffts;
        cfar_total(s,1) = cfar_total(s,1)+(half_N-2*margin)*(Nchirp-2*(Td+Gd));
        cfar_total(s,2) = cfar_total(s,2)+roi_cfar_cut_count;

        %% Map the sparse FFT output back to global range bins
        sparse_lookup = zeros(half_N,1);
        sparse_lookup(expanded_bins) = 1:numel(expanded_bins);

        %% One-to-one target matching and AoA estimation
        trial_hit = false(num_methods,Ntargets);
        trial_estimates = NaN(num_methods,Ntargets,num_metrics);

        for method = 1:num_methods
            if method==1
                peaks = dense_peaks;
            else
                peaks = sparse_peaks;
            end

            used_peaks = false(size(peaks,1),1);
            for k = 1:Ntargets
                [idx,used_peaks] = match_peak(peaks,used_peaks,true_range_bins(k),true_doppler_bins(k),data.RD_power);
                if idx==0
                    continue;
                end
                r = peaks(idx,1);
                d = peaks(idx,2);
                if method==1
                    rx_vector = squeeze(data.dense_fft(r,:,d));
                else
                    local_r = sparse_lookup(r);
                    assert(local_r>0);
                    rx_vector = squeeze(sparse_fft(local_r,:,d));
                end
                [azimuth,elevation] = planar_aoa(rx_vector,lambda,d_az,d_el,NangleFFT_az,NangleFFT_el);
                trial_hit(method,k) = true;
                trial_estimates(method,k,1) = range_axis(r);
                trial_estimates(method,k,2) = velocity_axis(d);
                trial_estimates(method,k,3) = azimuth;
                trial_estimates(method,k,4) = elevation;

                truth = [target_range(k),target_velocity(k),target_angle(k),target_elevation(k)];
                for j = 1:num_metrics
                    err = trial_estimates(method,k,j)-truth(j);
                    errors_all(s,method,j) = errors_all(s,method,j)+err^2;
                end
            end
        end

        for k = 1:Ntargets
            if trial_hit(2,k)
                assert(trial_hit(1,k),"Sparse hit not present in dense reference");
            end
            if all(trial_hit(:,k))
                common_hits(s) = common_hits(s)+1;
                for method = 1:num_methods
                    truth = [target_range(k),target_velocity(k),target_angle(k),target_elevation(k)];
                    for j = 1:num_metrics
                        err = trial_estimates(method,k,j)-truth(j);
                        errors_common(s,method,j) = errors_common(s,method,j)+err^2;
                    end
                end
                for j = 1:num_metrics
                    delta = abs(trial_estimates(1,k,j)-trial_estimates(2,k,j));
                    max_estimate_difference = max(max_estimate_difference,delta);
                    assert(delta<=1e-9,"Dense and Sparse common detections disagree");
                end
            end
        end
        hits(s,1) = hits(s,1)+sum(trial_hit(1,:));
        hits(s,2) = hits(s,2)+sum(trial_hit(2,:));
    end
    fprintf("Scene %d/%d done: %s | Dense hits %d | Sparse hits %d | Common %d\n",s,num_scenes,scene_names(s),hits(s,1),hits(s,2),common_hits(s));
end

%% Results: both all-detected and common-detection conditional RMSE
rows = num_scenes*num_methods;
Scene = strings(rows,1);
Method = strings(rows,1);
SNR_dB = zeros(rows,1);
ROIThreshold_dB = zeros(rows,1);
Trials = num_trials*ones(rows,1);
TargetEvents = num_trials*Ntargets*ones(rows,1);
MatchedTargets = zeros(rows,1);
MissedTargets = zeros(rows,1);
CommonMatchedTargets = zeros(rows,1);
MeanDopplerFFTs = zeros(rows,1);
MeanCFARCUTs = zeros(rows,1);
RangeRMSE_m = NaN(rows,1);
VelocityRMSE_mps = NaN(rows,1);
AzimuthRMSE_deg = NaN(rows,1);
ElevationRMSE_deg = NaN(rows,1);
CommonRangeRMSE_m = NaN(rows,1);
CommonVelocityRMSE_mps = NaN(rows,1);
CommonAzimuthRMSE_deg = NaN(rows,1);
CommonElevationRMSE_deg = NaN(rows,1);

row = 0;
for s = 1:num_scenes
    for method = 1:num_methods
        row = row+1;
        Scene(row) = scene_names(s);
        SNR_dB(row) = scene_snr(s);
        ROIThreshold_dB(row) = scene_threshold(s);
        if method==1
            Method(row) = "Dense";
        else
            Method(row) = "Balanced Sparse (actual FFT/CFAR)";
        end
        MatchedTargets(row) = hits(s,method);
        MissedTargets(row) = num_trials*Ntargets-hits(s,method);
        CommonMatchedTargets(row) = common_hits(s);
        MeanDopplerFFTs(row) = fft_total(s,method)/num_trials;
        MeanCFARCUTs(row) = cfar_total(s,method)/num_trials;
        if hits(s,method)>0
            RangeRMSE_m(row) = sqrt(errors_all(s,method,1)/hits(s,method));
            VelocityRMSE_mps(row) = sqrt(errors_all(s,method,2)/hits(s,method));
            AzimuthRMSE_deg(row) = sqrt(errors_all(s,method,3)/hits(s,method));
            ElevationRMSE_deg(row) = sqrt(errors_all(s,method,4)/hits(s,method));
        end
        if common_hits(s)>0
            CommonRangeRMSE_m(row) = sqrt(errors_common(s,method,1)/common_hits(s));
            CommonVelocityRMSE_mps(row) = sqrt(errors_common(s,method,2)/common_hits(s));
            CommonAzimuthRMSE_deg(row) = sqrt(errors_common(s,method,3)/common_hits(s));
            CommonElevationRMSE_deg(row) = sqrt(errors_common(s,method,4)/common_hits(s));
        end
    end
end

results = table(Scene,Method,SNR_dB,ROIThreshold_dB,Trials,TargetEvents,MatchedTargets,MissedTargets,CommonMatchedTargets,MeanDopplerFFTs,MeanCFARCUTs,RangeRMSE_m,VelocityRMSE_mps,AzimuthRMSE_deg,ElevationRMSE_deg,CommonRangeRMSE_m,CommonVelocityRMSE_mps,CommonAzimuthRMSE_deg,CommonElevationRMSE_deg);
assert(all(MatchedTargets+MissedTargets==TargetEvents));
assert(all(MatchedTargets(2:2:end)<=MatchedTargets(1:2:end)));
assert(all(CommonMatchedTargets(1:2:end)==CommonMatchedTargets(2:2:end)));

%% Save
output_dir = "output/phase2/16_dense_sparse_accuracy_validation";
if ~isfolder(output_dir)
    mkdir(output_dir);
end
csv_file = fullfile(output_dir,"dense_sparse_rmse_and_misses.csv");
txt_file = fullfile(output_dir,"dense_sparse_accuracy_validation.txt");
writetable(results,csv_file);

file = fopen(txt_file,"w");
assert(file~=-1,"Cannot write Group 5E TXT");
fprintf(file,"=== Group 5E: Actual Dense vs Sparse Accuracy Validation ===\n");
fprintf(file,"Independent trials per SNR scenario = %d\n",num_trials);
fprintf(file,"Target events per scenario = %d\n",num_trials*Ntargets);
fprintf(file,"Common-detection maximum estimate difference = %.3e\n",max_estimate_difference);
fprintf(file,"CFAR alpha = %.6f (original nominal, not calibrated 4.8)\n",Ntraining*(Pfa^(-1/Ntraining)-1));
fprintf(file,"Matching tolerance = +/-1 range bin and +/-1 Doppler bin, strongest unused 3x3 peak\n\n");
for i = 1:height(results)
    fprintf(file,"%s | %s | matched %d/%d | missed %d | common %d | Range %.4f m | Velocity %.4f m/s | Az %.4f deg | El %.4f deg\n",Scene(i),Method(i),MatchedTargets(i),TargetEvents(i),MissedTargets(i),CommonMatchedTargets(i),CommonRangeRMSE_m(i),CommonVelocityRMSE_mps(i),CommonAzimuthRMSE_deg(i),CommonElevationRMSE_deg(i));
end
fprintf(file,"\n[Limitations]\n");
fprintf(file,"Actual Sparse Doppler FFT and ROI CFAR run on the same simulated scene as Dense.\n");
fprintf(file,"Conditional RMSE is calculated only for matched targets and must be read together with misses.\n");
fprintf(file,"Ground truth is assigned using a +/-1-bin local peak rule; low-SNR noise hits are possible.\n");
fprintf(file,"One-to-one matching selects the strongest available peak within each target window.\n");
fprintf(file,"This test does not measure hardware latency, energy or overall Pd.\n");
fprintf(file,"No figures created. Group 5C remains separate mask-based statistical analysis.\n");
fprintf(file,"\nOverall code validation = PASS\n");
fclose(file);

fprintf("\n--- Group 5E Validation ---\n");
fprintf("Actual Sparse FFT / ROI CFAR / AoA = PASS\n");
fprintf("One-to-one matching / missing targets = PASS\n");
fprintf("Common-detection Dense vs Sparse estimate agreement = PASS\n");
fprintf("Group 5E code validation passed.\n");
fprintf("Results CSV: %s\n",csv_file);
fprintf("Summary: %s\n",txt_file);

%% Local helper: strongest unused grouped peak within +/-1 bin
function [index,used] = match_peak(peaks,used,true_r,true_d,power_map)
    index = 0;
    best_power = -inf;
    for i = 1:size(peaks,1)
        if used(i) || abs(peaks(i,1)-true_r)>1 || abs(peaks(i,2)-true_d)>1
            continue;
        end
        candidate_power = power_map(peaks(i,1),peaks(i,2));
        if candidate_power>best_power
            best_power = candidate_power;
            index = i;
        end
    end
    if index>0
        used(index) = true;
    end
end

%% Local helper: same 2x2 planar 2D angle FFT as src/angle/aoa_fft.m
function [azimuth,elevation] = planar_aoa(rx_vector,lambda,d_az,d_el,Naz,Nel)
    rx_vector = rx_vector(:);
    assert(numel(rx_vector)==4);
    rx_matrix = [rx_vector(1),rx_vector(2);rx_vector(3),rx_vector(4)];
    az_axis = (-Naz/2:Naz/2-1)*lambda/(Naz*d_az);
    el_axis = (-Nel/2:Nel/2-1)*lambda/(Nel*d_el);
    [az_grid,el_grid] = meshgrid(az_axis,el_axis);
    power = abs(fftshift(fft2(rx_matrix,Nel,Naz))).^2;
    power(az_grid.^2+el_grid.^2>1) = 0;
    [~,peak_idx] = max(power(:));
    [el_idx,az_idx] = ind2sub(size(power),peak_idx);
    elevation = asind(el_axis(el_idx));
    az_argument = az_axis(az_idx)/cosd(elevation);
    azimuth = asind(max(-1,min(1,az_argument)));
end
