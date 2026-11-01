
function [selected_bins, blocks, expanded_bins, dropped_bins] = roi_budget_selection(candidate_bins, range_energy, half_N, margin, max_rois, max_fft_count, Nrx)

%% Input Validation

candidate_bins = unique(candidate_bins(:));
range_energy = range_energy(:);

assert(length(range_energy) == half_N);
assert(all(isfinite(range_energy)));
assert(all(candidate_bins >= 1 & candidate_bins <= half_N));
assert(max_rois >= 0 && max_rois == floor(max_rois));
assert(max_fft_count >= 0 && max_fft_count == floor(max_fft_count));
assert(Nrx > 0);


%% Rank Candidates by Range Energy

if isempty(candidate_bins)

    selected_bins = zeros(0,1);
    dropped_bins = zeros(0,1);

    [blocks, expanded_bins] = roi_merge_blocks(selected_bins, half_N, margin);

    return;

end

scores = range_energy(candidate_bins);

ranking = [-scores, candidate_bins];

[~, order] = sortrows(ranking, [1 2]);

ranked_bins = candidate_bins(order);


%% Greedy Selection With Hard FFT Budget

selected_bins = zeros(0,1);

expanded_mask = false(half_N,1);

for i = 1:length(ranked_bins)

    if length(selected_bins) >= max_rois
        break;
    end

    r = ranked_bins(i);

    first_bin = max(1, r-margin);
    last_bin = min(half_N, r+margin);

    new_mask = expanded_mask;

    new_mask(first_bin:last_bin) = true;

    required_ffts = nnz(new_mask) * Nrx;

    if required_ffts <= max_fft_count

        selected_bins(end+1,1) = r;

        expanded_mask = new_mask;

    end

end


%% Merge Selected ROI Into Contiguous Blocks

selected_bins = sort(selected_bins);

[blocks, expanded_bins] = roi_merge_blocks(selected_bins, half_N, margin);

dropped_bins = setdiff(candidate_bins, selected_bins);


%% Final Hard-Cap Validation

actual_ffts = length(expanded_bins) * Nrx;

assert(length(selected_bins) <= max_rois);
assert(actual_ffts <= max_fft_count);
assert(isequal(expanded_bins, find(expanded_mask)));
assert(all(ismember(selected_bins, candidate_bins)));
assert(isempty(intersect(selected_bins, dropped_bins)));

end
