function print_roi_cfar_margin_selection(result, half_N, Nrx)

%% ROI CFAR Margin 선택 결과 출력

fprintf("\n--- ROI CFAR Margin Selection ---\n");

% 선택된 ROI 및 CFAR Margin 정보
fprintf("Original ROI bins = %d\n", length(result.roi_bins));
fprintf("CFAR range margin = %d bins per side\n", result.roi_range_margin);

% Margin 확장 후 처리할 Range Bin 정보
fprintf("ROI + margin bins = %d / %d\n", result.roi_processing_count, half_N);
fprintf("Additional margin bins = %d\n", result.roi_extra_margin_bins);

% 필요한 Doppler FFT 실행 횟수
fprintf("Required Doppler FFTs = %d\n", result.roi_processing_count * Nrx);

end