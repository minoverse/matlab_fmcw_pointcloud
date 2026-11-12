function result = run_radar_pipeline()
%% 모듈화된 FMCW Radar 파이프라인
root = fileparts(mfilename("fullpath"));
cd(root);
addpath(genpath(fullfile(root,"src")));
cfg = radar_config();
snr_db = 20;
[adc_cube, t_fast] = simulate_adc_cube_process(cfg, snr_db);
plot_first_chirp_iq(adc_cube, t_fast);
[range_data, range_axis, range_data_no_window] = range_fft_process(adc_cube, cfg);
plot_range_fft(range_data, range_data_no_window, range_axis);
range_data_for_doppler = static_clutter_process(range_data);
apply_clutter_removal = false;
[range_doppler_rx, RD_power, velocity_axis] = doppler_fft_process(range_data_for_doppler, cfg, apply_clutter_removal);
plot_doppler_fft(RD_power, range_axis, velocity_axis, cfg);
[cfar_map, noise_map, det_range_idx, det_doppler_idx] = cfar2d_process(RD_power, cfg.cfar);
plot_cfar_detections(RD_power, det_range_idx, det_doppler_idx, range_axis, velocity_axis);
[peak_range_idx, peak_doppler_idx] = peak_grouping_process(RD_power, det_range_idx, det_doppler_idx);
[estimated_azimuth, estimated_elevation, estimated_snr] = aoa_fft_process(range_doppler_rx, RD_power, noise_map, peak_range_idx, peak_doppler_idx, cfg);
plot_aoa_detections(estimated_azimuth, estimated_elevation);
pointcloud = build_pointcloud_process(peak_range_idx, peak_doppler_idx, range_axis, velocity_axis, estimated_azimuth, estimated_elevation, estimated_snr);
validation = validate_ground_truth_process(pointcloud,cfg);
print_ground_truth_results(validation);
save_pointcloud_outputs(pointcloud);
result = struct("config",cfg,"pointcloud",pointcloud,"RD_power",RD_power,"cfar_map",cfar_map,"validation",validation);
end
