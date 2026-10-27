%% Peak grouping

peak_range_idx = [];

peak_doppler_idx = [];

peak_power = [];

for i = 1:length(det_range_idx)

    r = det_range_idx(i);

    d = det_doppler_idx(i);

    if r <= 1 || r >= size(RD_power, 1)
        continue;
    end

    if d <= 1 || d >= size(RD_power, 2)
        continue;
    end

    local_area = RD_power(r-1:r+1, d-1:d+1);

    CUT = RD_power(r, d);

    if CUT == max(local_area(:))

        peak_range_idx(end+1) = r;

        peak_doppler_idx(end+1) = d;

        peak_power(end+1) = CUT;

    end

end

fprintf("\n--- Peak Grouping ---\n");

fprintf("CFAR cells = %d\n", length(det_range_idx));

fprintf("Grouped peaks = %d\n", length(peak_range_idx));

for i = 1:length(peak_range_idx)

    peak_range = range_axis(peak_range_idx(i));

    peak_velocity = velocity_axis(peak_doppler_idx(i));

    fprintf("Peak %d: Range = %.2f m, Velocity = %.2f m/s\n", i, peak_range, peak_velocity);

end