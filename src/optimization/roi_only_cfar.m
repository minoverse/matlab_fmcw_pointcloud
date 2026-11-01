
%% ROI-Only CA-CFAR

%% CFAR Parameters

range_half_width = Tr + Gr;

doppler_half_width = Td + Gd;

roi_total_cells = (2 * range_half_width + 1) * (2 * doppler_half_width + 1);

roi_guard_cells = (2 * Gr + 1) * (2 * Gd + 1);

roi_training_cells = roi_total_cells - roi_guard_cells;

roi_alpha = roi_training_cells * (Pfa^(-1 / roi_training_cells) - 1);


%% Training Mask

training_mask_roi = true(2 * range_half_width + 1, 2 * doppler_half_width + 1);

center_r = range_half_width + 1;

center_d = doppler_half_width + 1;

training_mask_roi(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = false;


%% Map Original Range Indices to Sparse Rows

range_lookup = zeros(half_N, 1);

range_lookup(active_range_bins) = 1:length(active_range_bins);


%% ROI CFAR Results

roi_cfar_local = false(length(roi_bins), Nchirp);

roi_noise_local = zeros(length(roi_bins), Nchirp);

roi_cfar_cut_count = 0;


%% CFAR on ROI CUTs Only

for i = 1:length(roi_bins)

    r = roi_bins(i);

    if r <= range_half_width || r > half_N - range_half_width
        continue;
    end

    neighborhood = r-range_half_width:r+range_half_width;

    local_rows = range_lookup(neighborhood);

    assert(all(local_rows > 0), "Missing CFAR training margin");

    for d = doppler_half_width+1:Nchirp-doppler_half_width

        window = sparse_RD_power(local_rows, d-doppler_half_width:d+doppler_half_width);

        training_cells_roi = window(training_mask_roi);

        noise_power_roi = mean(training_cells_roi);

        roi_noise_local(i, d) = noise_power_roi;

        threshold_roi = roi_alpha * noise_power_roi;

        CUT = sparse_RD_power(range_lookup(r), d);

        if CUT > threshold_roi
            roi_cfar_local(i, d) = true;
        end

        roi_cfar_cut_count = roi_cfar_cut_count + 1;

    end

end


%% Detection Indices in Original Range Coordinates

[roi_local_det_idx, roi_det_doppler_idx] = find(roi_cfar_local);

roi_det_range_idx = roi_bins(roi_local_det_idx);


%% Peak Grouping Using Sparse RD Power

roi_peak_range_idx = [];

roi_peak_doppler_idx = [];

for i = 1:length(roi_det_range_idx)

    r = roi_det_range_idx(i);

    d = roi_det_doppler_idx(i);

    if r <= 1 || r >= half_N
        continue;
    end

    if d <= 1 || d >= Nchirp
        continue;
    end

    local_rows = range_lookup(r-1:r+1);

    assert(all(local_rows > 0), "Missing peak grouping margin");

    local_area = sparse_RD_power(local_rows, d-1:d+1);

    CUT = sparse_RD_power(range_lookup(r), d);

    if CUT == max(local_area(:))

        roi_peak_range_idx(end+1) = r;

        roi_peak_doppler_idx(end+1) = d;

    end

end


%% Validation

valid_roi_bins = roi_bins(roi_bins > range_half_width & roi_bins <= half_N-range_half_width);

expected_cuts = length(valid_roi_bins) * (Nchirp - 2 * doppler_half_width);

assert(roi_training_cells == Ntraining);

assert(roi_cfar_cut_count == expected_cuts);


%% Results

fprintf("\n--- ROI-Only CFAR ---\n");

fprintf("Training cells = %d\n", roi_training_cells);

fprintf("ROI CFAR CUTs = %d\n", roi_cfar_cut_count);

fprintf("ROI CFAR training accesses = %d\n", roi_cfar_cut_count * roi_training_cells);

fprintf("ROI CFAR detected cells = %d\n", length(roi_det_range_idx));

fprintf("ROI grouped peaks = %d\n", length(roi_peak_range_idx));
