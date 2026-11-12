function print_roi_only_cfar(roi_stats)

%% ROI-Only CFAR 결과 출력

fprintf("\n--- ROI-Only CFAR ---\n");

% Training Cell 개수
fprintf("Training cells = %d\n", roi_stats.training_cells);

% 실제 CFAR CUT 처리량
fprintf("ROI CFAR CUTs = %d\n", roi_stats.cfar_cuts);

% Training Cell 접근 횟수
fprintf("ROI CFAR training accesses = %d\n", roi_stats.training_accesses);

% CFAR 검출 및 Peak Grouping 결과
fprintf("ROI CFAR detected cells = %d\n", roi_stats.detected_cells);
fprintf("ROI grouped peaks = %d\n", roi_stats.grouped_peaks);

end