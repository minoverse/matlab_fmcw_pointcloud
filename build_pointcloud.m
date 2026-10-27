%% Build Point Cloud

Npoints = length(peak_range_idx);

detected_range = zeros(Npoints, 1);

detected_velocity = zeros(Npoints, 1);

x = zeros(Npoints, 1);

y = zeros(Npoints, 1);

for i = 1:Npoints

    detected_range(i) = range_axis(peak_range_idx(i));

    detected_velocity(i) = velocity_axis(peak_doppler_idx(i));

    x(i) = detected_range(i) * cosd(estimated_angles(i));

    y(i) = detected_range(i) * sind(estimated_angles(i));

end

pointcloud = table(x, y, detected_range, detected_velocity, estimated_angles, estimated_snr);

pointcloud.Properties.VariableNames = {'x', 'y', 'range', 'radial_velocity', 'angle', 'snr'};
fprintf("\n--- Point Cloud ---\n");

disp(pointcloud);

%% CSV
if ~exist("output", "dir")
    mkdir("output");
end

writetable(pointcloud, "output/pointcloud.csv");

%% Plot
figure;

scatter(x, y, 70, "filled");

xlabel("Forward x [m]");

ylabel("Lateral y [m]");

title("Radar Point Cloud");

grid on;

axis equal;

if ~exist("plots/06_pointcloud", "dir")
    mkdir("plots/06_pointcloud");
end

exportgraphics(gcf, "plots/06_pointcloud/pointcloud.png", "Resolution", 300);