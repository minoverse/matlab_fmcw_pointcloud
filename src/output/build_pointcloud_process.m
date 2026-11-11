
function pointcloud = build_pointcloud_process( ...
    peak_range_idx, peak_doppler_idx, ...
    range_axis, velocity_axis, ...
    estimated_azimuth, estimated_elevation, estimated_snr)

%% Point Cloud 배열 초기화

Npoints = length(peak_range_idx);

detected_range = zeros(Npoints, 1);
detected_velocity = zeros(Npoints, 1);

x = zeros(Npoints, 1);
y = zeros(Npoints, 1);
z = zeros(Npoints, 1);


%% Range , Velocity , 3D 좌표 계산

for i = 1:Npoints

    detected_range(i) = range_axis(peak_range_idx(i));

    detected_velocity(i) = velocity_axis(peak_doppler_idx(i));


    range_value = detected_range(i);

    azimuth = estimated_azimuth(i);
    elevation = estimated_elevation(i);


    % 구면 좌표 직교 좌표 변환
    x(i) = range_value * cosd(elevation) * cosd(azimuth);

    y(i) = range_value * cosd(elevation) * sind(azimuth);

    z(i) = range_value * sind(elevation);

end


%% Point Cloud Table 생성

pointcloud = table( x, y, z, detected_range, detected_velocity, ...
    estimated_azimuth, estimated_elevation, estimated_snr);

pointcloud.Properties.VariableNames = {'x', 'y', 'z', 'range', 'radial_velocity', ...
    'azimuth', 'elevation', 'snr'};

end
