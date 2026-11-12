function [stress_rows,stress_candidate_counts,stress_fft_caps] = roi_budget_regression(half_N,Nrx,Tr,Gr)


%% Group 4B - Synthetic Workload Stress Test

stress_candidate_counts = [0, 1, 5, 10, 25, 50, 100, 250, 500];

stress_fft_caps = [0, 80, 100, 200, 400];

stress_top_k = 6;

stress_margin = Tr + Gr;

stress_energy = (1:half_N)' + 10*sin((1:half_N)');

stress_rows = zeros(length(stress_candidate_counts)*length(stress_fft_caps), 6);

row = 0;


%% Test Different ROI Densities

for i = 1:length(stress_candidate_counts)

    candidate_count = stress_candidate_counts(i);

    if candidate_count == 0

        candidates = zeros(0,1);

    else

        candidates = unique(round(linspace(1, half_N, candidate_count)))';

    end


    for j = 1:length(stress_fft_caps)

        fft_cap = stress_fft_caps(j);


        %% Apply Bounded ROI Selection

        [kept, test_blocks, test_expanded, dropped] = roi_budget_selection(candidates, stress_energy, half_N, stress_margin, stress_top_k, fft_cap, Nrx);


        %% Verify Hard Limits

        fft_count = length(test_expanded) * Nrx;

        assert(length(kept) <= stress_top_k);

        assert(fft_count <= fft_cap);

        assert(all(ismember(kept, candidates)));

        assert(isempty(intersect(kept, dropped)));

        assert(isequal(sort([kept; dropped]), candidates));


        %% Verify CFAR Margin Coverage

        expected_mask = false(half_N,1);

        for k = 1:length(kept)

            r = kept(k);

            first_bin = max(1, r-stress_margin);
            last_bin = min(half_N, r+stress_margin);

            expected_mask(first_bin:last_bin) = true;

        end

        assert(isequal(test_expanded, find(expected_mask)));


        %% Verify Block Length

        total_block_bins = sum(test_blocks(:,2) - test_blocks(:,1) + 1);

        assert(total_block_bins == length(test_expanded));


        %% Save Counters

        row = row + 1;

        stress_rows(row,:) = [candidate_count, fft_cap, length(kept), size(test_blocks,1), length(test_expanded), fft_count];

    end

end


%% Equal-Energy Ranking Test

equal_energy = ones(half_N,1);

[tie_bins, ~, ~, ~] = roi_budget_selection([10;20], equal_energy, half_N, stress_margin, 1, 100, Nrx);

assert(isequal(tie_bins,10));


%% Zero-Budget Test

[zero_bins, zero_blocks, zero_expanded, ~] = roi_budget_selection((1:half_N)', stress_energy, half_N, stress_margin, 10, 0, Nrx);

assert(isempty(zero_bins));
assert(isempty(zero_blocks));
assert(isempty(zero_expanded));


%% Full-Budget Test

[full_bins, full_blocks, full_expanded, full_dropped] = roi_budget_selection((1:half_N)', stress_energy, half_N, stress_margin, half_N, half_N*Nrx, Nrx);

assert(length(full_bins) == half_N);
assert(isequal(full_blocks, [1,half_N]));
assert(length(full_expanded) == half_N);
assert(isempty(full_dropped));


%% Final Result

fprintf("\n--- Group 4B Synthetic Stress Test ---\n");

fprintf("Candidate densities tested = %d\n", length(stress_candidate_counts));

fprintf("FFT budgets tested = %d\n", length(stress_fft_caps));

fprintf("Total density/budget cases = %d\n", size(stress_rows,1));

fprintf("Top-K limit validation = PASS\n");
fprintf("Hard FFT budget validation = PASS\n");
fprintf("Margin coverage validation = PASS\n");
fprintf("Equal-energy ranking = PASS\n");
fprintf("Zero-budget handling = PASS\n");
fprintf("Full-budget handling = PASS\n");
fprintf("Synthetic workload stress test = PASS\n");


end
