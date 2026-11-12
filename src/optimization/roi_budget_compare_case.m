function result = roi_budget_compare_case(method, candidate_bins, data, radar, top_k, fft_cap)

    %% Radar 파라미터

    half_N = radar.half_N;
    Nrx = radar.Nrx;
    Nchirp = radar.Nchirp;

    Tr = radar.Tr;
    Gr = radar.Gr;

    margin = Tr + Gr;

    candidate_bins = unique(candidate_bins(:));


    %% ROI 선택 방식에 따른 처리

    if method == "Bins"

        [selected_bins, processing_blocks, expanded_bins, dropped_bins] = roi_budget_selection(candidate_bins,data.range_energy,half_N,margin,top_k,fft_cap,Nrx);

        [source_blocks, ~] = roi_merge_blocks(candidate_bins,half_N,0);

        [selected_source_blocks, ~] = roi_merge_blocks(selected_bins,half_N,0);

        selected_units = length(selected_bins);

    elseif method == "Blocks"

        [selected_bins, selected_source_blocks, processing_blocks, expanded_bins, dropped_bins, source_blocks] = roi_block_budget_selection(candidate_bins,data.range_energy,half_N,margin,top_k,fft_cap,Nrx);

        selected_units = size(selected_source_blocks,1);

    else

        error("Unknown ROI selection method");

    end


    %% 선택된 Block에 Doppler FFT 수행

    [block_doppler_rx, block_RD_power, block_fft_count] = block_doppler_processing(data.range_data,processing_blocks,false);


    %% FFT Hard Cap 검증

    assert(selected_units <= top_k);

    assert(block_fft_count <= fft_cap);

    assert(block_fft_count == length(expanded_bins)*Nrx);

    assert(length(selected_bins)+length(dropped_bins) == length(candidate_bins));


    %% Dense FFT와 Block FFT 결과 비교

    if ~isempty(expanded_bins)

        dense_selected_fft = data.dense_fft(expanded_bins,:,:);

        fft_error = max(abs(block_doppler_rx(:)-dense_selected_fft(:)));

        fft_tolerance = 1e-10*max(1,max(abs(dense_selected_fft(:))));

        assert(fft_error <= fft_tolerance);


        dense_selected_power = data.RD_power(expanded_bins,:);

        power_error = max(abs(block_RD_power(:)-dense_selected_power(:)));

        power_tolerance = 1e-10*max(1,max(abs(dense_selected_power(:))));

        assert(power_error <= power_tolerance);

    end


    %% 선택된 ROI에만 CA-CFAR 수행

    roi_bins = selected_bins;

    active_range_bins = expanded_bins;

    sparse_RD_power = block_RD_power;

    [roi_cfar_local, roi_noise_local, roi_peak_range_idx, roi_peak_doppler_idx, roi_cfar_cut_count, roi_stats] = roi_only_cfar_process(roi_bins, active_range_bins, sparse_RD_power, radar);

    print_roi_only_cfar(roi_stats);


    %% Dense CFAR와 Sparse CFAR 결과 비교

    expected_cfar = data.cfar_map(roi_bins,:);

    assert(isequal(roi_cfar_local,expected_cfar));

    if ~isempty(roi_bins)

        expected_noise = data.noise_map(roi_bins,:);

        noise_error = max(abs(roi_noise_local(:)-expected_noise(:)));

        noise_tolerance = 1e-10*max(1,max(abs(expected_noise(:))));

        assert(noise_error <= noise_tolerance);

    end


    %% Dense Peak와 Sparse Peak 결과 비교

    [dense_local_r,dense_local_d] = find(expected_cfar);

    expected_peaks = zeros(0,2);

    for i = 1:length(dense_local_r)

        r = roi_bins(dense_local_r(i));
        d = dense_local_d(i);

        if r <= 1 || r >= half_N || d <= 1 || d >= Nchirp
            continue;
        end

        local_area = data.RD_power(r-1:r+1,d-1:d+1);

        if data.RD_power(r,d) == max(local_area(:))

            expected_peaks(end+1,:) = [r,d];

        end

    end

    sparse_peaks = [roi_peak_range_idx(:),roi_peak_doppler_idx(:)];

    assert(isequal(sortrows(expected_peaks),sortrows(sparse_peaks)));


    %% 실제 Target 검출 결과 비교

    bounded_map = false(half_N,Nchirp);

    bounded_map(roi_bins,:) = roi_cfar_local;

    num_targets = length(data.target_ranges);

    dense_hits = zeros(num_targets,1);
    bounded_hits = zeros(num_targets,1);

    for k = 1:num_targets

        [~,r0] = min(abs(data.range_axis-data.target_ranges(k)));

        [~,d0] = min(abs(data.velocity_axis-data.target_velocities(k)));

        r_near = max(1,r0-1):min(half_N,r0+1);

        d_near = max(1,d0-1):min(Nchirp,d0+1);

        dense_hits(k) = any(any(data.cfar_map(r_near,d_near)));

        bounded_hits(k) = any(any(bounded_map(r_near,d_near)));

        assert(bounded_hits(k) <= dense_hits(k));

    end


    %% ROI Budget 비교 결과 반환

    result.method = method;

    result.candidate_bins = length(candidate_bins);
    result.candidate_blocks = size(source_blocks,1);

    result.selected_units = selected_units;
    result.selected_bins = length(selected_bins);

    result.selected_source_blocks = size(selected_source_blocks,1);
    result.processing_blocks = size(processing_blocks,1);

    result.dropped_bins = length(dropped_bins);

    result.expanded_bins = length(expanded_bins);
    result.fft_count = block_fft_count;

    result.cfar_cuts = roi_cfar_cut_count;
    result.cfar_detections = nnz(roi_cfar_local);

    result.dense_hits = dense_hits;
    result.bounded_hits = bounded_hits;

    result.fft_reduction = 100*(1-block_fft_count/(half_N*Nrx));

end