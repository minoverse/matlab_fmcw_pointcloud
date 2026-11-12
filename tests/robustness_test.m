function result = robustness_test(test_case, near_zero_notch)
%% FMCW Robustness Test - 모듈 함수 기반

if nargin < 1
    test_case = 6;
end
if nargin < 2
    near_zero_notch = true;
end

root = fileparts(fileparts(mfilename("fullpath")));
cd(root);
addpath(genpath(fullfile(root, "src")));

cfg = radar_config();
apply_clutter_removal = false;

if test_case == 1
    cfg.target_amplitude = [1.0, 0.2];
end

if test_case == 2
    cfg.target_range = [30, 30.45];
    cfg.target_velocity = [10, 10];
    cfg.target_angle = [15, -20];
    cfg.target_elevation = [10, -8];
    cfg.target_amplitude = [1.0, 0.7];
end

if test_case == 3
    cfg.target_range = [30, 55, 20];
    cfg.target_velocity = [10, -5, 0];
    cfg.target_angle = [15, -20, 0];
    cfg.target_elevation = [10, -8, 0];
    cfg.target_amplitude = [1.0, 0.7, 2.0];
    apply_clutter_removal = true;
end

if test_case == 5
    cfg.target_range = [30, 55, 40];
    cfg.target_velocity = [10, -5, 0.6];
    cfg.target_angle = [15, -20, 0];
    cfg.target_elevation = [10, -8, 5];
    cfg.target_amplitude = [1.0, 0.7, 1.0];
    apply_clutter_removal = near_zero_notch;
end

if test_case == 6
    cfg.target_elevation = [20, -15];
end

% Target 파라미터 변경 후 파생값 재계산
cfg.Ntargets = length(cfg.target_range);
cfg.tau = 2 * cfg.target_range / cfg.c;
cfg.fb = cfg.S .* cfg.tau;
cfg.fd = 2 * cfg.target_velocity / cfg.lambda;

[adc_cube,t_fast] = simulate_adc_cube_process(cfg,20);
plot_first_chirp_iq(adc_cube,t_fast);

if test_case == 4
    rx_gain = [1.0, 0.9, 1.1, 0.85];
    rx_phase = [0, 8, -5, 12];
    for rx = 1:cfg.Nrx
        phase_error = exp(1j * rx_phase(rx) * pi / 180);
        adc_cube(:,rx,:) = adc_cube(:,rx,:) * rx_gain(rx) * phase_error;
    end
end

[range_data,range_axis,range_data_no_window] = range_fft_process(adc_cube,cfg);
plot_range_fft(range_data,range_data_no_window,range_axis);
range_data_for_doppler = static_clutter_process(range_data);
[range_doppler_rx,RD_power,velocity_axis] = doppler_fft_process(range_data_for_doppler,cfg,apply_clutter_removal);
plot_doppler_fft(RD_power,range_axis,velocity_axis,cfg);
[cfar_map,noise_map,det_range_idx,det_doppler_idx] = cfar2d_process(RD_power,cfg.cfar);
plot_cfar_detections(RD_power,det_range_idx,det_doppler_idx,range_axis,velocity_axis);
[peak_range_idx,peak_doppler_idx] = peak_grouping_process(RD_power,det_range_idx,det_doppler_idx);
[estimated_azimuth,estimated_elevation,estimated_snr] = aoa_fft_process(range_doppler_rx,RD_power,noise_map,peak_range_idx,peak_doppler_idx,cfg);
plot_aoa_detections(estimated_azimuth,estimated_elevation);
pointcloud = build_pointcloud_process(peak_range_idx,peak_doppler_idx,range_axis,velocity_axis,estimated_azimuth,estimated_elevation,estimated_snr);
validation = validate_ground_truth_process(pointcloud,cfg);
print_ground_truth_results(validation);

% 원래 Robustness 결과 파일명 유지
if ~isfolder("plots/robustness")
    mkdir("plots/robustness");
end
if ~isfolder("output")
    mkdir("output");
end

if test_case == 5
    if near_zero_notch
        filename = "plots/robustness/test_5_notch_on_pointcloud.png";
        csv_name = "output/robustness_test_5_notch_on.csv";
    else
        filename = "plots/robustness/test_5_notch_off_pointcloud.png";
        csv_name = "output/robustness_test_5_notch_off.csv";
    end
else
    filename = "plots/robustness/test_" + test_case + "_pointcloud.png";
    csv_name = "output/robustness_test_" + test_case + ".csv";
end

figure;
scatter3(pointcloud.x,pointcloud.y,pointcloud.z,80,"filled");
xlabel("Forward x [m]");
ylabel("Lateral y [m]");
zlabel("Height z [m]");
title("Robustness Test Point Cloud");
grid on;
axis equal;
view(3);
exportgraphics(gcf,filename,"Resolution",300);
writetable(pointcloud,csv_name);

result = struct("pointcloud",pointcloud,"validation",validation,"config",cfg);
end
