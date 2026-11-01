
%% Lightweight Pre-Doppler ROI Gating

%% Positive Range FFT

range_positive = range_data(1:half_N, :, :);


%% Noncoherent Range Energy

range_energy = squeeze(sum(sum(abs(range_positive).^2, 3), 2));


%% Median-Based Threshold

roi_threshold_db = 6;

roi_scale = 10^(roi_threshold_db / 10);

roi_noise_floor = median(range_energy);

roi_threshold = roi_scale * roi_noise_floor;


%% ROI Candidate Selection

roi_mask = range_energy > roi_threshold;

roi_bins = find(roi_mask);

roi_count = length(roi_bins);

roi_fraction = roi_count / half_N;


%% Target Validation

target_roi_selected = false(Ntargets, 1);

for k = 1:Ntargets

    [~, target_bin] = min(abs(range_axis - target_range(k)));

    target_roi_selected(k) = roi_mask(target_bin);

end


%% Assertions

assert(length(range_energy) == half_N);
assert(length(roi_mask) == half_N);

assert(roi_count > 0, "No ROI bins selected");
assert(roi_count < half_N, "ROI did not reduce range bins");

assert(all(target_roi_selected), "At least one target missed by ROI");


%% Results

fprintf("\n--- Pre-Doppler ROI Gating ---\n");

fprintf("Threshold = %.2f dB above median\n", roi_threshold_db);

fprintf("Median energy = %.4f\n", roi_noise_floor);

fprintf("Selected range bins = %d / %d\n", roi_count, half_N);

fprintf("Selected fraction = %.2f%%\n", roi_fraction * 100);

for k = 1:Ntargets

    fprintf("Target %d ROI selected = %d\n", k, target_roi_selected(k));

end


%% Plot

figure;

plot(range_axis, range_energy);

hold on;

yline(roi_threshold, "r--");

plot(range_axis(roi_bins), range_energy(roi_bins), "go");

xlabel("Range [m]");

ylabel("Noncoherent Range Energy");

title("Pre-Doppler ROI Candidate Selection");

legend("Range Energy", "Threshold", "Selected ROI");

grid on;

%% Save ROI Selection Plot

output_dir = "output/phase2/01_baseline_roi_gating";

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

plot_name = sprintf("pre_doppler_roi_energy_threshold_%gdB.png", roi_threshold_db);

exportgraphics(gcf, fullfile(output_dir, plot_name), "Resolution", 300);

fprintf("ROI plot saved: %s\n", fullfile(output_dir, plot_name));
