
%% Baseline Cost Breakdown

%% Input Validation

assert(size(adc_cube, 1) == Ns);
assert(size(adc_cube, 2) == Nrx);
assert(size(adc_cube, 3) == Nchirp);

assert(size(range_data, 1) == Ns);
assert(size(range_doppler_rx, 1) == Ns);

assert(size(RD_power, 1) == half_N);
assert(size(RD_power, 2) == Nchirp);

assert(roi_count == length(roi_bins));


%% Range FFT

range_fft_count = Nrx * Nchirp;

range_fft_work = range_fft_count * Ns * log2(Ns);

% Extra FFT for window comparison in range_fft.m
range_fft_diagnostic_count = range_fft_count;


%% Doppler FFT

% Current implementation processes all range bins
doppler_fft_current = Ns * Nrx;

% Fair Dense processes positive range bins only
doppler_fft_fair = half_N * Nrx;

% Projected Sparse processes ROI bins only
doppler_fft_sparse = roi_count * Nrx;

doppler_work_per_fft = Nchirp * log2(Nchirp);

doppler_work_current = doppler_fft_current * doppler_work_per_fft;

doppler_work_fair = doppler_fft_fair * doppler_work_per_fft;

doppler_work_sparse = doppler_fft_sparse * doppler_work_per_fft;

projected_roi_reduction = 100 * (1 - doppler_work_sparse / doppler_work_fair);


%% CFAR

range_width = 2 * (Tr + Gr) + 1;

doppler_width = 2 * (Td + Gd) + 1;

total_cells_check = range_width * doppler_width;

guard_cells_check = (2 * Gr + 1) * (2 * Gd + 1);

training_cells_check = total_cells_check - guard_cells_check;

valid_range_cuts = half_N - 2 * (Tr + Gr);

valid_doppler_cuts = Nchirp - 2 * (Td + Gd);

cfar_cut_count = valid_range_cuts * valid_doppler_cuts;

cfar_training_accesses = cfar_cut_count * Ntraining;


%% AoA

aoa_fft_count = Npeaks;

aoa_grid_points = NangleFFT_az * NangleFFT_el;

aoa_work = aoa_fft_count * aoa_grid_points * (log2(NangleFFT_az) + log2(NangleFFT_el));


%% Memory Size Proxies

complex_bytes = 16;
real_bytes = 8;

adc_cube_bytes = numel(adc_cube) * complex_bytes;

range_cube_bytes = numel(range_data) * complex_bytes;

doppler_cube_current_bytes = numel(range_doppler_rx) * complex_bytes;

doppler_cube_fair_bytes = half_N * Nrx * Nchirp * complex_bytes;

doppler_cube_sparse_bytes = roi_count * Nrx * Nchirp * complex_bytes;

rd_power_bytes = numel(RD_power) * real_bytes;


%% Logical FFT I/O Proxies

range_fft_io_bytes = 2 * range_fft_count * Ns * complex_bytes;

doppler_io_current_bytes = 2 * doppler_fft_current * Nchirp * complex_bytes;

doppler_io_fair_bytes = 2 * doppler_fft_fair * Nchirp * complex_bytes;

doppler_io_sparse_bytes = 2 * doppler_fft_sparse * Nchirp * complex_bytes;


%% ROI Gating Work Proxy

roi_energy_samples = half_N * Nrx * Nchirp;


%% Counter Validation

assert(range_fft_count == 256);

assert(doppler_fft_current == 4000);
assert(doppler_fft_fair == 2000);
assert(doppler_fft_sparse == roi_count * Nrx);

assert(range_width == 21);
assert(doppler_width == 11);

assert(total_cells == 231);
assert(guard_cells == 15);
assert(Ntraining == 216);

assert(total_cells_check == total_cells);
assert(guard_cells_check == guard_cells);
assert(training_cells_check == Ntraining);

assert(valid_range_cuts == 480);
assert(valid_doppler_cuts == 54);

assert(cfar_cut_count == 25920);
assert(cfar_training_accesses == 5598720);

assert(aoa_fft_count == length(peak_range_idx));


%% Results

fprintf("\n--- Baseline Cost Breakdown ---\n");

fprintf("\nRange FFT:\n");
fprintf("DSP Range FFT count = %d\n", range_fft_count);
fprintf("Diagnostic Range FFT count = %d\n", range_fft_diagnostic_count);
fprintf("Range FFT work proxy = %.0f\n", range_fft_work);

fprintf("\nDoppler FFT:\n");
fprintf("Current Dense FFT count = %d\n", doppler_fft_current);
fprintf("Fair Dense FFT count = %d\n", doppler_fft_fair);
fprintf("Projected Sparse FFT count = %d\n", doppler_fft_sparse);

fprintf("Current Dense work proxy = %.0f\n", doppler_work_current);
fprintf("Fair Dense work proxy = %.0f\n", doppler_work_fair);
fprintf("Projected Sparse work proxy = %.0f\n", doppler_work_sparse);

fprintf("Projected ROI Doppler work reduction vs Fair Dense = %.2f%%\n", projected_roi_reduction);

fprintf("\nCFAR:\n");
fprintf("Total cells = %d\n", total_cells);
fprintf("Guard cells = %d\n", guard_cells);
fprintf("Training cells = %d\n", Ntraining);
fprintf("CFAR CUT count = %d\n", cfar_cut_count);
fprintf("Training accesses = %d\n", cfar_training_accesses);

fprintf("\nAoA:\n");
fprintf("AoA FFT count = %d\n", aoa_fft_count);
fprintf("AoA work proxy = %.0f\n", aoa_work);

fprintf("\nMemory Size Proxies:\n");
fprintf("ADC cube = %.2f MiB\n", adc_cube_bytes / 2^20);
fprintf("Range cube = %.2f MiB\n", range_cube_bytes / 2^20);
fprintf("Current Doppler cube = %.2f MiB\n", doppler_cube_current_bytes / 2^20);
fprintf("Fair Dense Doppler cube = %.2f MiB\n", doppler_cube_fair_bytes / 2^20);
fprintf("Projected Sparse Doppler cube = %.2f MiB\n", doppler_cube_sparse_bytes / 2^20);
fprintf("RD power map = %.2f MiB\n", rd_power_bytes / 2^20);

fprintf("\nLogical FFT I/O Proxies:\n");
fprintf("Range FFT I/O = %.2f MiB\n", range_fft_io_bytes / 2^20);
fprintf("Current Doppler I/O = %.2f MiB\n", doppler_io_current_bytes / 2^20);
fprintf("Fair Dense Doppler I/O = %.2f MiB\n", doppler_io_fair_bytes / 2^20);
fprintf("Projected Sparse Doppler I/O = %.2f MiB\n", doppler_io_sparse_bytes / 2^20);

fprintf("\nROI Energy Samples = %d\n", roi_energy_samples);

fprintf("\nBaseline counter validation = PASS\n");
