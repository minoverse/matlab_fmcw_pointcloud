%% 2D CA-CFAR

%% CFAR parameters
Tr = 8;
Td = 4;

Gr = 2;
Gd = 1;

Pfa = 1e-5;

%% Number of training cells
total_cells = (2*(Tr + Gr) + 1) * (2*(Td + Gd) + 1);

guard_cells = (2*Gr + 1) * (2*Gd + 1);

Ntraining = total_cells - guard_cells;

%% Threshold factor
alpha = Ntraining * (Pfa^(-1/Ntraining) - 1);

%% Maps
cfar_map = false(size(RD_power));

noise_map = zeros(size(RD_power));

[Nrange, Ndoppler] = size(RD_power);

%% Move CUT through Range-Doppler map
for r = Tr + Gr + 1 : Nrange - (Tr + Gr)

    for d = Td + Gd + 1 : Ndoppler - (Td + Gd)

        window = RD_power(r-(Tr+Gr):r+(Tr+Gr), d-(Td+Gd):d+(Td+Gd));

        training_mask = true(size(window));

        center_r = Tr + Gr + 1;

        center_d = Td + Gd + 1;

        training_mask(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = false;

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

%% Detected cells
[det_range_idx, det_doppler_idx] = find(cfar_map);

detected_ranges = range_axis(det_range_idx);

detected_velocities = velocity_axis(det_doppler_idx);

%% Print
fprintf("\n--- 2D CA-CFAR ---\n");

fprintf("Training cells = %d\n", Ntraining);

fprintf("Pfa = %.1e\n", Pfa);

fprintf("Detected cells = %d\n", length(det_range_idx));

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

plot(detected_velocities, detected_ranges, "ro", "MarkerSize", 7, "LineWidth", 1.5);

hold off;

%% Save
if ~exist("plots/04_cfar", "dir")
    mkdir("plots/04_cfar");
end

exportgraphics(gcf, "plots/04_cfar/cfar_detections.png", "Resolution", 300);