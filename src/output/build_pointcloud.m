%% Build 3D Point Cloud

Npoints = length(peak_range_idx);

detected_range = zeros(Npoints, 1);

detected_velocity = zeros(Npoints, 1);

x = zeros(Npoints, 1);

y = zeros(Npoints, 1);

z = zeros(Npoints, 1);


for i = 1:Npoints

    detected_range(i) = range_axis(peak_range_idx(i));

    detected_velocity(i) = velocity_axis(peak_doppler_idx(i));


    range_value = detected_range(i);

    azimuth = estimated_azimuth(i);

    elevation = estimated_elevation(i);


    x(i) = range_value * cosd(elevation) * cosd(azimuth);

    y(i) = range_value * cosd(elevation) * sind(azimuth);

    z(i) = range_value * sind(elevation);

end


%% Point cloud table

pointcloud = table(x, y, z, detected_range, detected_velocity, estimated_azimuth, estimated_elevation, estimated_snr);


pointcloud.Properties.VariableNames = {'x', 'y', 'z', 'range', 'radial_velocity', 'azimuth', 'elevation', 'snr'};


fprintf("\n--- 3D Point Cloud ---\n");

disp(pointcloud);


%% CSV

if ~exist("output", "dir")

    mkdir("output");

end


writetable(pointcloud, "output/pointcloud.csv");


%% Plot

figure;

scatter3(x, y, z, 80, "filled");

xlabel("Forward x [m]");

ylabel("Lateral y [m]");

zlabel("Height z [m]");

title("Radar 3D Point Cloud");

grid on;

axis equal;

view(3);


%% Save

if ~exist("plots/06_pointcloud", "dir")

    mkdir("plots/06_pointcloud");

end


exportgraphics(gcf, "plots/06_pointcloud/pointcloud.png", "Resolution", 300);