
%% ROI CFAR Margin Selection

%% Check ROI

roi_bins = unique(roi_bins(:));

assert(~isempty(roi_bins), "No ROI bins selected");
assert(all(roi_bins >= 1));
assert(all(roi_bins <= half_N));


%% CFAR Range Margin

roi_range_margin = Tr + Gr;

roi_margin_mask = false(half_N, 1);


%% Expand ROI With Margin

for i = 1:length(roi_bins)

    r = roi_bins(i);

    first_bin = max(1, r - roi_range_margin);

    last_bin = min(half_N, r + roi_range_margin);

    roi_margin_mask(first_bin:last_bin) = true;

end


%% Range Bins Required for Doppler FFT

active_range_bins = find(roi_margin_mask);

roi_processing_count = length(active_range_bins);

roi_extra_margin_bins = roi_processing_count - length(roi_bins);


%% Validation

assert(all(roi_margin_mask(roi_bins)));

assert(roi_processing_count >= length(roi_bins));


%% Results

fprintf("\n--- ROI CFAR Margin Selection ---\n");

fprintf("Original ROI bins = %d\n", length(roi_bins));

fprintf("CFAR range margin = %d bins per side\n", roi_range_margin);

fprintf("ROI + margin bins = %d / %d\n", roi_processing_count, half_N);

fprintf("Additional margin bins = %d\n", roi_extra_margin_bins);

fprintf("Required Doppler FFTs = %d\n", roi_processing_count * Nrx);

