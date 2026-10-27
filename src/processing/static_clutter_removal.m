%% Static / DC clutter handling

if ~exist("apply_clutter_removal", "var")
    apply_clutter_removal = false;
end

range_data_for_doppler = range_data;