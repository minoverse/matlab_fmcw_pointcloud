function print_roi_candidate_selection(roi, half_N)

%% Results

fprintf("\n--- Pre-Doppler ROI Gating ---\n");

fprintf("Threshold = %.2f dB above median\n", roi.roi_threshold_db);
fprintf("Median energy = %.4f\n", roi.roi_noise_floor);

fprintf("Selected range bins = %d / %d\n", roi.roi_count, half_N);
fprintf("Selected fraction = %.2f%%\n", roi.roi_fraction * 100);

for k = 1:length(roi.target_roi_selected)

    fprintf("Target %d ROI selected = %d\n", k, roi.target_roi_selected(k));

end

end