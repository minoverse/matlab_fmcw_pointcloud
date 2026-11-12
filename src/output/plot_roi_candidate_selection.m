function plot_roi_candidate_selection(roi, range_axis, output_dir)

%% Plot

figure;

plot(range_axis, roi.range_energy);

hold on;

yline(roi.roi_threshold, "r--");

plot(range_axis(roi.roi_bins), roi.range_energy(roi.roi_bins), "go");

xlabel("Range [m]");
ylabel("Noncoherent Range Energy");

title("Pre-Doppler ROI Candidate Selection");

legend("Range Energy", "Threshold", "Selected ROI");

grid on;


%% Save Plot

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

plot_name = sprintf("pre_doppler_roi_energy_threshold_%gdB.png", roi.roi_threshold_db);

exportgraphics(gcf, fullfile(output_dir, plot_name), "Resolution", 300);

fprintf("ROI plot saved: %s\n", fullfile(output_dir, plot_name));

end