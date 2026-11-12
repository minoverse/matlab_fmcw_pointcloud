function roi = roi_candidate_selection_process(range_data, half_N, range_axis, target_range, roi_threshold_db)

%% Positive Range FFT

range_positive = range_data(1:half_N, :, :);


%% Noncoherent Range Energy

range_energy = squeeze(sum(sum(abs(range_positive).^2, 3), 2));


%% Median-Based Threshold

roi_scale = 10^(roi_threshold_db / 10);

roi_noise_floor = median(range_energy);

roi_threshold = roi_scale * roi_noise_floor;


%% ROI Candidate Selection

roi_mask = range_energy > roi_threshold;

roi_bins = find(roi_mask);

roi_count = length(roi_bins);

roi_fraction = roi_count / half_N;


%% Target Validation

Ntargets = length(target_range);

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


%% Return Results

roi.range_energy = range_energy;

roi.roi_threshold_db = roi_threshold_db;
roi.roi_noise_floor = roi_noise_floor;
roi.roi_threshold = roi_threshold;

roi.roi_mask = roi_mask;
roi.roi_bins = roi_bins;

roi.roi_count = roi_count;
roi.roi_fraction = roi_fraction;

roi.target_roi_selected = target_roi_selected;

end