
function cost = baseline_cost_process(cfg, data)

    %% 입력 검증

    Ns = cfg.Ns;
    Nrx = cfg.Nrx;
    Nchirp = cfg.Nchirp;
    half_N = floor(Ns / 2);

    roi_count = data.roi_count;

    assert(size(data.adc_cube, 1) == Ns);
    assert(size(data.adc_cube, 2) == Nrx);
    assert(size(data.adc_cube, 3) == Nchirp);

    assert(size(data.range_data, 1) == Ns);
    assert(size(data.range_doppler_rx, 1) == Ns);

    assert(size(data.RD_power, 1) == half_N);
    assert(size(data.RD_power, 2) == Nchirp);

    assert(roi_count == length(data.roi_bins));
    assert(data.Npeaks == length(data.peak_range_idx));


    %% Range FFT 연산량

    cost.range_fft_count = Nrx * Nchirp;
    cost.range_fft_work = cost.range_fft_count * Ns * log2(Ns);

    % Window 비교용 추가 FFT
    cost.range_fft_diagnostic_count = cost.range_fft_count;


    %% Doppler FFT 연산량

    % 현재 Dense: 전체 Range Bin
    cost.doppler_fft_current = Ns * Nrx;

    % Fair Dense: 양의 Range Bin
    cost.doppler_fft_fair = half_N * Nrx;

    % 예상 Sparse: ROI Bin
    cost.doppler_fft_sparse = roi_count * Nrx;

    cost.doppler_work_per_fft = Nchirp * log2(Nchirp);

    cost.doppler_work_current = cost.doppler_fft_current * cost.doppler_work_per_fft;
    cost.doppler_work_fair = cost.doppler_fft_fair * cost.doppler_work_per_fft;
    cost.doppler_work_sparse = cost.doppler_fft_sparse * cost.doppler_work_per_fft;

    cost.projected_roi_reduction = 100 * (1 - cost.doppler_work_sparse / cost.doppler_work_fair);


    %% CFAR 연산량

    Tr = cfg.cfar.Tr;
    Td = cfg.cfar.Td;
    Gr = cfg.cfar.Gr;
    Gd = cfg.cfar.Gd;

    cost.range_width = 2 * (Tr + Gr) + 1;
    cost.doppler_width = 2 * (Td + Gd) + 1;

    cost.total_cells = cost.range_width * cost.doppler_width;
    cost.guard_cells = (2 * Gr + 1) * (2 * Gd + 1);
    cost.Ntraining = cost.total_cells - cost.guard_cells;

    cost.valid_range_cuts = half_N - 2 * (Tr + Gr);
    cost.valid_doppler_cuts = Nchirp - 2 * (Td + Gd);

    cost.cfar_cut_count = cost.valid_range_cuts * cost.valid_doppler_cuts;
    cost.cfar_training_accesses = cost.cfar_cut_count * cost.Ntraining;


    %% AoA 연산량

    cost.aoa_fft_count = data.Npeaks;

    cost.aoa_grid_points = cfg.NangleFFT_az * cfg.NangleFFT_el;

    cost.aoa_work = cost.aoa_fft_count * cost.aoa_grid_points * (log2(cfg.NangleFFT_az) + log2(cfg.NangleFFT_el));


    %% 메모리 크기 추정

    complex_bytes = 16;
    real_bytes = 8;

    cost.adc_cube_bytes = numel(data.adc_cube) * complex_bytes;
    cost.range_cube_bytes = numel(data.range_data) * complex_bytes;
    cost.doppler_cube_current_bytes = numel(data.range_doppler_rx) * complex_bytes;

    cost.doppler_cube_fair_bytes = half_N * Nrx * Nchirp * complex_bytes;
    cost.doppler_cube_sparse_bytes = roi_count * Nrx * Nchirp * complex_bytes;

    cost.rd_power_bytes = numel(data.RD_power) * real_bytes;


    %% FFT 논리 I/O 추정

    cost.range_fft_io_bytes = 2 * cost.range_fft_count * Ns * complex_bytes;

    cost.doppler_io_current_bytes = 2 * cost.doppler_fft_current * Nchirp * complex_bytes;
    cost.doppler_io_fair_bytes = 2 * cost.doppler_fft_fair * Nchirp * complex_bytes;
    cost.doppler_io_sparse_bytes = 2 * cost.doppler_fft_sparse * Nchirp * complex_bytes;


    %% ROI 검사량

    cost.roi_energy_samples = half_N * Nrx * Nchirp;


    %% 계산 일관성 확인

    assert(cost.total_cells == cost.guard_cells + cost.Ntraining);
    assert(cost.doppler_fft_sparse == roi_count * Nrx);
    assert(cost.cfar_training_accesses == cost.cfar_cut_count * cost.Ntraining);

end
