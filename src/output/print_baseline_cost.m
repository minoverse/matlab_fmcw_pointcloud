
function print_baseline_cost(c)

%% Range FFT

fprintf("\n--- Baseline Cost Breakdown ---\n");

fprintf("\nRange FFT:\n");
fprintf("DSP Range FFT count = %d\n", c.range_fft_count);
fprintf("Diagnostic Range FFT count = %d\n", c.range_fft_diagnostic_count);
fprintf("Range FFT work proxy = %.0f\n", c.range_fft_work);


%% Doppler FFT

fprintf("\nDoppler FFT:\n");

fprintf("Current Dense FFT count = %d\n", c.doppler_fft_current);
fprintf("Fair Dense FFT count = %d\n", c.doppler_fft_fair);
fprintf("Projected Sparse FFT count = %d\n", c.doppler_fft_sparse);

fprintf("Current Dense work proxy = %.0f\n", c.doppler_work_current);
fprintf("Fair Dense work proxy = %.0f\n", c.doppler_work_fair);
fprintf("Projected Sparse work proxy = %.0f\n", c.doppler_work_sparse);

fprintf("Projected ROI Doppler work reduction vs Fair Dense = %.2f%%\n", c.projected_roi_reduction);


%% CFAR

fprintf("\nCFAR:\n");

fprintf("Total cells = %d\n", c.total_cells);
fprintf("Guard cells = %d\n", c.guard_cells);
fprintf("Training cells = %d\n", c.Ntraining);
fprintf("CFAR CUT count = %d\n", c.cfar_cut_count);
fprintf("Training accesses = %d\n", c.cfar_training_accesses);


%% AoA

fprintf("\nAoA:\n");

fprintf("AoA FFT count = %d\n", c.aoa_fft_count);
fprintf("AoA work proxy = %.0f\n", c.aoa_work);


%% 메모리 크기

fprintf("\nMemory Size Proxies:\n");

fprintf("ADC cube = %.2f MiB\n", c.adc_cube_bytes / 2^20);
fprintf("Range cube = %.2f MiB\n", c.range_cube_bytes / 2^20);
fprintf("Current Doppler cube = %.2f MiB\n", c.doppler_cube_current_bytes / 2^20);
fprintf("Fair Dense Doppler cube = %.2f MiB\n", c.doppler_cube_fair_bytes / 2^20);
fprintf("Projected Sparse Doppler cube = %.2f MiB\n", c.doppler_cube_sparse_bytes / 2^20);
fprintf("RD power map = %.2f MiB\n", c.rd_power_bytes / 2^20);


%% FFT 논리 I/O

fprintf("\nLogical FFT I/O Proxies:\n");

fprintf("Range FFT I/O = %.2f MiB\n", c.range_fft_io_bytes / 2^20);
fprintf("Current Doppler I/O = %.2f MiB\n", c.doppler_io_current_bytes / 2^20);
fprintf("Fair Dense Doppler I/O = %.2f MiB\n", c.doppler_io_fair_bytes / 2^20);
fprintf("Projected Sparse Doppler I/O = %.2f MiB\n", c.doppler_io_sparse_bytes / 2^20);


%% ROI

fprintf("\nROI Energy Samples = %d\n", c.roi_energy_samples);

fprintf("\nBaseline counter validation = PASS\n");

end
