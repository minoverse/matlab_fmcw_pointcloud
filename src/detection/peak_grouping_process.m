
function [peak_range_idx, peak_doppler_idx, peak_power] = peak_grouping_process(RD_power, det_range_idx, det_doppler_idx)

%% Initialize peaks
peak_range_idx = [];
peak_doppler_idx = [];
peak_power = [];


%% Local peak grouping
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

end
