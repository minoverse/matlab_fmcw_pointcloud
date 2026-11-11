
function save_pointcloud_outputs(pointcloud)

%% Point Cloud 결과 출력

fprintf("\n--- 3D Point Cloud ---\n");

disp(pointcloud);


%% CSV 저장

if ~exist("output", "dir")
    mkdir("output");
end

writetable(pointcloud, "output/pointcloud.csv");


%% 3D Point Cloud 시각화

figure;

scatter3( pointcloud.x, pointcloud.y, pointcloud.z, 80, "filled");

xlabel("Forward x [m]");
ylabel("Lateral y [m]");
zlabel("Height z [m]");

title("Radar 3D Point Cloud");

grid on;
axis equal;
view(3);


%% PNG 저장

if ~exist("plots/06_pointcloud", "dir")
    mkdir("plots/06_pointcloud");
end

exportgraphics(gcf, ...
    "plots/06_pointcloud/pointcloud.png", "Resolution", 300);

end
