function print_selective_doppler_processing(active_range_bins, half_N, sparse_fft_count, sparse_fft_work)

%% 선택적 Doppler 처리 결과 출력

fprintf("\n--- Selective Doppler Processing ---\n");

% 처리할 Range Bin 개수
fprintf("Active range bins = %d / %d\n", length(active_range_bins), half_N);

% 실제 Doppler FFT 실행 횟수
fprintf("Doppler FFTs executed = %d\n", sparse_fft_count);

% Doppler FFT 연산량 추정치
fprintf("Doppler FFT work proxy = %.0f\n", sparse_fft_work);

end