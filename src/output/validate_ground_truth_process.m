
function result = validate_ground_truth_process(pointcloud, cfg)

%% 결과 배열 초기화

Npoints = height(pointcloud);

result.gt_idx = zeros(Npoints, 1);
result.range_error = zeros(Npoints, 1);
result.velocity_error = zeros(Npoints, 1);
result.azimuth_error = zeros(Npoints, 1);
result.elevation_error = zeros(Npoints, 1);


%% Ground Truth와 검출 결과 비교

for i = 1:Npoints

    % Range 기준 가장 가까운 Target 선택
    [~, gt_idx] = min(abs(cfg.target_range - pointcloud.range(i)));

    result.gt_idx(i) = gt_idx;

    % 측정값과 실제값의 오차 계산
    result.range_error(i) = pointcloud.range(i) - cfg.target_range(gt_idx);

    result.velocity_error(i) = pointcloud.radial_velocity(i) - cfg.target_velocity(gt_idx);

    result.azimuth_error(i) = pointcloud.azimuth(i) - cfg.target_angle(gt_idx);

    result.elevation_error(i) = pointcloud.elevation(i) - cfg.target_elevation(gt_idx);

end


%% 이론적 Radar 성능 계산

result.range_resolution = cfg.range_resolution;

result.velocity_resolution = cfg.lambda / (2 * cfg.Nchirp * cfg.Tchirp);

result.max_velocity = cfg.lambda / (4 * cfg.Tchirp);

end
