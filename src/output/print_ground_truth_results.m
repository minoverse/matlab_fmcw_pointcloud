
function print_ground_truth_results(result)

%% Ground Truth 오차 출력

fprintf("\n--- Ground Truth Validation ---\n");

for i = 1:length(result.gt_idx)

    fprintf("\nDetection %d -> Target %d\n", i, result.gt_idx(i));

    fprintf("Range error = %.3f m\n", result.range_error(i));

    fprintf("Velocity error = %.3f m/s\n", result.velocity_error(i));

    fprintf("Azimuth error = %.3f deg\n", result.azimuth_error(i));

    fprintf("Elevation error = %.3f deg\n", result.elevation_error(i));

end


%% 이론적 성능 출력

fprintf("\n--- Theoretical Validation ---\n");

fprintf("Range resolution = %.3f m\n", result.range_resolution);

fprintf("Velocity resolution = %.3f m/s\n", result.velocity_resolution);

fprintf("Max unambiguous velocity = %.3f m/s\n", result.max_velocity);

end
