
function [roi_mask, roi_bins, expanded_bins] = roi_threshold_gate(range_positive, threshold_db, Tr, Gr)

%% Range Energy

range_energy = squeeze(sum(sum(abs(range_positive).^2, 3), 2));

half_N = size(range_positive, 1);


%% Median Threshold

noise_floor = median(range_energy);

threshold = noise_floor * 10^(threshold_db / 10);


%% ROI Selection

roi_mask = range_energy > threshold;

roi_bins = find(roi_mask);


%% CFAR Range Margin

margin = Tr + Gr;

expanded_mask = false(half_N, 1);

for i = 1:length(roi_bins)

    r = roi_bins(i);

    first_bin = max(1, r - margin);

    last_bin = min(half_N, r + margin);

    expanded_mask(first_bin:last_bin) = true;

end

expanded_bins = find(expanded_mask);

end

