
function [blocks, expanded_bins] = roi_merge_blocks(roi_bins, half_N, margin)

%% Check Input

roi_bins = unique(roi_bins(:));

assert(half_N > 0);
assert(margin >= 0);
assert(all(roi_bins >= 1));
assert(all(roi_bins <= half_N));
assert(all(roi_bins == floor(roi_bins)));


%% Empty ROI

blocks = zeros(0, 2);
expanded_bins = zeros(0, 1);

if isempty(roi_bins)
    return;
end


%% Expand and Merge

blocks = zeros(length(roi_bins), 2);

block_count = 0;

for i = 1:length(roi_bins)

    r = roi_bins(i);

    first_bin = max(1, r - margin);
    last_bin = min(half_N, r + margin);

    if block_count == 0

        block_count = 1;
        blocks(block_count, :) = [first_bin, last_bin];

    elseif first_bin <= blocks(block_count, 2) + 1

        blocks(block_count, 2) = max(blocks(block_count, 2), last_bin);

    else

        block_count = block_count + 1;
        blocks(block_count, :) = [first_bin, last_bin];

    end

end

blocks = blocks(1:block_count, :);


%% Build Unique Processing Bin List

total_bins = sum(blocks(:, 2) - blocks(:, 1) + 1);

expanded_bins = zeros(total_bins, 1);

position = 1;

for i = 1:block_count

    first_bin = blocks(i, 1);
    last_bin = blocks(i, 2);

    count = last_bin - first_bin + 1;

    expanded_bins(position:position+count-1) = (first_bin:last_bin)';

    position = position + count;

end


%% Validation

assert(length(expanded_bins) == total_bins);
assert(all(diff(expanded_bins) > 0));
assert(all(ismember(roi_bins, expanded_bins)));

end

