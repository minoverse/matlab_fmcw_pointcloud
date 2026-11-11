
function plot_aoa_detections(estimated_azimuth, estimated_elevation)

%% AoA 검출 결과 시각화

figure;

scatter(estimated_azimuth, estimated_elevation, 80, "filled");

xlabel("Azimuth [deg]");
ylabel("Elevation [deg]");

title("Azimuth / Elevation Detections");

grid on;

xlim([-90 90]);
ylim([-60 60]);


%% 결과 이미지 저장

output_dir = "plots/05_aoa";

if ~isfolder(output_dir)
    mkdir(output_dir);
end

exportgraphics(gcf, fullfile(output_dir, "aoa_detections.png"), ...
    "Resolution", 300);

end
