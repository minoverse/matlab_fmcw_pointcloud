
function [cfar_map, noise_map, det_range_idx, det_doppler_idx, stats] = ...
    cfar2d_process(RD_power, cfg)

%% CFAR parameters
Tr = cfg.Tr;
Td = cfg.Td;

Gr = cfg.Gr;
Gd = cfg.Gd;

Pfa = cfg.Pfa;


%% Number of training cells
total_cells = (2*(Tr+Gr)+1) * (2*(Td+Gd)+1);

guard_cells = (2*Gr+1) * (2*Gd+1);

Ntraining = total_cells - guard_cells;


%% Threshold factor
alpha = Ntraining * (Pfa^(-1/Ntraining) - 1);


%% Initialize maps
cfar_map = false(size(RD_power));

noise_map = zeros(size(RD_power));

[Nrange, Ndoppler] = size(RD_power);


%% 2D CA-CFAR
for r = Tr+Gr+1 : Nrange-(Tr+Gr)

    for d = Td+Gd+1 : Ndoppler-(Td+Gd)

        window = RD_power( ...
            r-(Tr+Gr):r+(Tr+Gr), ...
            d-(Td+Gd):d+(Td+Gd));

        training_mask = true(size(window));

        center_r = Tr+Gr+1;
        center_d = Td+Gd+1;

        training_mask(center_r-Gr:center_r+Gr, ...
            center_d-Gd:center_d+Gd) = false;

        training_cells = window(training_mask);

        noise_power_cfar = mean(training_cells);

        noise_map(r, d) = noise_power_cfar;

        threshold = alpha * noise_power_cfar;

        CUT = RD_power(r, d);

        if CUT > threshold
            cfar_map(r, d) = true;
        end

    end

end


%% Detection indices
[det_range_idx, det_doppler_idx] = find(cfar_map);


%% Processing statistics
stats.total_cells = total_cells;
stats.guard_cells = guard_cells;
stats.Ntraining = Ntraining;
stats.alpha = alpha;
stats.Pfa = Pfa;
stats.detected_cells = length(det_range_idx);

end
