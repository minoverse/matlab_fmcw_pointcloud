
function plot_cfar_detections(RD_power, det_range_idx, det_doppler_idx, range_axis, velocity_axis)

%% Detection coordinates
detected_ranges = range_axis(det_range_idx);

detected_velocities = velocity_axis(det_doppler_idx);


%% Range-Doppler magnitude
RD_mag = sqrt(RD_power);

RD_db = 20*log10(RD_mag / max(RD_mag(:)) + eps);


%% Plot
figure;

imagesc(velocity_axis, range_axis, RD_db);

axis xy;

xlabel("Radial Velocity [m/s]");
ylabel("Range [m]");

title("2D CA-CFAR Detections");

colorbar;

ylim([0 70]);
caxis([-60 0]);

hold on;

plot(detected_velocities, detected_ranges, ...
    "ro", "MarkerSize", 7, "LineWidth", 1.5);

hold off;


%% Save
output_dir = "plots/04_cfar";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

exportgraphics(gcf, fullfile(output_dir, "cfar_detections.png"), ...
    "Resolution", 300);

end
