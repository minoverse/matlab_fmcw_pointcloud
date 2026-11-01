
%% Group 4B-2 - Block Selection Regression

test_energy = zeros(half_N,1);

test_energy(100:102) = [100;95;90];

test_energy(300:301) = [70;60];

test_candidates = [100;101;102;300;301];

test_top_k = 2;

test_fft_cap = 180;


%% Existing Top-K Bins

[bin_selected, ~, bin_expanded, ~] = roi_budget_selection(test_candidates,test_energy,half_N,Tr+Gr,test_top_k,test_fft_cap,Nrx);


%% New Top-K Blocks

[block_selected, selected_source_blocks, processing_blocks, block_expanded, ~, source_blocks] = roi_block_budget_selection(test_candidates,test_energy,half_N,Tr+Gr,test_top_k,test_fft_cap,Nrx);


%% Compare Selection Behavior

assert(length(bin_selected) <= test_top_k);

assert(length(bin_expanded)*Nrx <= test_fft_cap);

assert(~any(ismember(bin_selected,[300;301])));

assert(any(ismember(block_selected,[300;301])));

assert(size(source_blocks,1) == 2);

assert(size(selected_source_blocks,1) == 2);

assert(size(processing_blocks,1) == 2);

assert(length(block_expanded)*Nrx == 180);


%% Equal-Energy Tie

equal_energy = ones(half_N,1);

[tie_selected, tie_blocks, ~, ~, ~, ~] = roi_block_budget_selection([20;21;100;101],equal_energy,half_N,Tr+Gr,1,100,Nrx);

assert(isequal(tie_selected,[20;21]));

assert(isequal(tie_blocks,[20,21]));


%% Zero FFT Budget

[zero_selected, zero_blocks, zero_processing, zero_expanded, ~, ~] = roi_block_budget_selection(test_candidates,test_energy,half_N,Tr+Gr,2,0,Nrx);

assert(isempty(zero_selected));
assert(isempty(zero_blocks));
assert(isempty(zero_processing));
assert(isempty(zero_expanded));


%% Zero Top-K Budget

[kzero_selected, ~, ~, kzero_expanded, ~, ~] = roi_block_budget_selection(test_candidates,test_energy,half_N,Tr+Gr,0,200,Nrx);

assert(isempty(kzero_selected));
assert(isempty(kzero_expanded));


%% Full Budget

[full_selected, full_selected_blocks, full_processing, full_expanded, full_dropped, ~] = roi_block_budget_selection((1:half_N)',equal_energy,half_N,Tr+Gr,half_N,half_N*Nrx,Nrx);

assert(length(full_selected) == half_N);
assert(size(full_selected_blocks,1) == 1);
assert(isequal(full_processing,[1,half_N]));
assert(length(full_expanded) == half_N);
assert(isempty(full_dropped));


%% Oversized Atomic Block

[oversize_selected, ~, ~, oversize_expanded, ~, ~] = roi_block_budget_selection((1:half_N)',equal_energy,half_N,Tr+Gr,6,200,Nrx);

assert(isempty(oversize_selected));
assert(isempty(oversize_expanded));


%% Results

fprintf("\n--- Group 4B-2 Synthetic Regression ---\n");

fprintf("Top-K Bins strong-target concentration = PASS\n");
fprintf("Top-K Blocks second-target coverage = PASS\n");
fprintf("Hard FFT budget = PASS\n");
fprintf("Equal-energy ranking = PASS\n");
fprintf("Zero-budget handling = PASS\n");
fprintf("Full-budget handling = PASS\n");
fprintf("Oversized block handling = PASS\n");
fprintf("Group 4B-2 synthetic regression = PASS\n");

