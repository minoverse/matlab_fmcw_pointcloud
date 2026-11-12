
function [selected_bins, selected_source_blocks, processing_blocks, expanded_bins, dropped_bins, source_blocks] = roi_block_budget_selection(candidate_bins, range_energy, half_N, margin, max_blocks, max_fft_count, Nrx)

    %% 입력값 검증

    candidate_bins = unique(candidate_bins(:));
    range_energy = range_energy(:);

    assert(length(range_energy) == half_N);
    assert(all(isfinite(range_energy)));
    assert(all(candidate_bins >= 1 & candidate_bins <= half_N));
    assert(max_blocks >= 0 && max_blocks == floor(max_blocks));
    assert(max_fft_count >= 0 && max_fft_count == floor(max_fft_count));
    assert(Nrx > 0);


    %% 인접한 후보 Bin을 Source Block으로 병합

    [source_blocks, ~] = roi_merge_blocks(candidate_bins, half_N, 0);

    num_blocks = size(source_blocks, 1);

    block_scores = zeros(num_blocks, 1);

    for i = 1:num_blocks

        first_bin = source_blocks(i, 1);
        last_bin = source_blocks(i, 2);

        block_scores(i) = max(range_energy(first_bin:last_bin));

    end


    %% 최대 에너지 기준으로 Block 우선순위 정렬

    ranking = [-block_scores, source_blocks(:,1)];

    [~, order] = sortrows(ranking, [1 2]);


    %% 연산량 제한을 적용한 Block 선택

    selected_block_mask = false(num_blocks, 1);

    selected_bin_mask = false(half_N, 1);

    processing_mask = false(half_N, 1);

    selected_count = 0;

    for i = 1:num_blocks

        if selected_count >= max_blocks
            break;
        end

        block_id = order(i);

        first_bin = source_blocks(block_id, 1);
        last_bin = source_blocks(block_id, 2);

        first_margin = max(1, first_bin-margin);
        last_margin = min(half_N, last_bin+margin);

        new_processing_mask = processing_mask;

        new_processing_mask(first_margin:last_margin) = true;

        required_ffts = nnz(new_processing_mask) * Nrx;

        if required_ffts <= max_fft_count

            selected_block_mask(block_id) = true;

            selected_bin_mask(first_bin:last_bin) = true;

            processing_mask = new_processing_mask;

            selected_count = selected_count + 1;

        end

    end


    %% 선택된 ROI 및 Processing Block 구성

    selected_bins = find(selected_bin_mask);

    selected_source_blocks = source_blocks(selected_block_mask, :);

    [processing_blocks, expanded_bins] = roi_merge_blocks(selected_bins, half_N, margin);

    dropped_bins = setdiff(candidate_bins, selected_bins);


    %% 연산량 제한 및 선택 결과 검증

    assert(size(selected_source_blocks,1) <= max_blocks);

    assert(length(expanded_bins)*Nrx <= max_fft_count);

    assert(isequal(expanded_bins, find(processing_mask)));

    assert(all(ismember(selected_bins, candidate_bins)));

    assert(isempty(intersect(selected_bins, dropped_bins)));

    assert(length(selected_bins)+length(dropped_bins) == length(candidate_bins));

end
