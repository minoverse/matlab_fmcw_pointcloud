function result = roi_cfar_margin_selection_process(roi_bins, half_N, Tr, Gr)

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


%% Return Results

result.roi_bins = roi_bins;
result.roi_range_margin = roi_range_margin;
result.roi_margin_mask = roi_margin_mask;

result.active_range_bins = active_range_bins;
result.roi_processing_count = roi_processing_count;
result.roi_extra_margin_bins = roi_extra_margin_bins;

end