%% Ground Truth Validation

fprintf("\n--- Ground Truth Validation ---\n");


for i = 1:Npoints

    [~, gt_idx] = min(abs(target_range - detected_range(i)));


    range_error = detected_range(i) - target_range(gt_idx);

    velocity_error = detected_velocity(i) - target_velocity(gt_idx);

    azimuth_error = estimated_azimuth(i) - target_angle(gt_idx);

    elevation_error = estimated_elevation(i) - target_elevation(gt_idx);


    fprintf("\nDetection %d -> Target %d\n", i, gt_idx);

    fprintf("Range error = %.3f m\n", range_error);

    fprintf("Velocity error = %.3f m/s\n", velocity_error);

    fprintf("Azimuth error = %.3f deg\n", azimuth_error);

    fprintf("Elevation error = %.3f deg\n", elevation_error);

end


%% Theoretical validation

velocity_resolution = lambda / (2 * Nchirp * Tchirp);

max_velocity = lambda / (4 * Tchirp);


fprintf("\n--- Theoretical Validation ---\n");

fprintf("Range resolution = %.3f m\n", range_resolution);

fprintf("Velocity resolution = %.3f m/s\n", velocity_resolution);

fprintf("Max unambiguous velocity = %.3f m/s\n", max_velocity);