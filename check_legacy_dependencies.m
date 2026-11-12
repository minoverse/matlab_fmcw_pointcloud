function [unused, blockers] = check_legacy_dependencies()
%% 삭제 후보 및 MATLAB 소스 내부 Script 의존성 점검
root = fileparts(mfilename("fullpath"));
legacy = [
    "config.m"
    "src/angle/aoa_fft.m"
    "src/detection/cfar2d.m"
    "src/detection/peak_grouping.m"
    "src/optimization/baseline_cost_breakdown.m"
    "src/optimization/roi_candidate_selection.m"
    "src/optimization/roi_cfar_margin_selection.m"
    "src/optimization/roi_only_cfar.m"
    "src/optimization/selective_doppler_processing.m"
    "src/output/build_pointcloud.m"
    "src/output/validate_ground_truth.m"
    "src/processing/doppler_fft.m"
    "src/processing/range_fft.m"
    "src/processing/static_clutter_removal.m"
    "src/simulation/simulate_adc_cube.m"
];
blockers = strings(0,1);
source_files = [dir(fullfile(root,"*.m")); dir(fullfile(root,"src","**","*.m")); dir(fullfile(root,"tests","**","*.m")); dir(fullfile(root,"experiments","**","*.m"))];
for i = 1:numel(source_files)
    f = fullfile(source_files(i).folder,source_files(i).name);
    if string(f) == string(fullfile(root,"check_legacy_dependencies.m"))
        continue;
    end
    src = fileread(f);
    src = regexprep(src,'(?m)^\s*%[^\n]*','');
    for j = 1:numel(legacy)
        [~,base] = fileparts(legacy(j));
        % 함수/파일 호출 이름까지 확인. 검사 대상 파일 자신은 제외.
        target = fullfile(root,legacy(j));
        if string(f)==string(target)
            continue;
        end
        pat = ['(?<![\w])' regexptranslate('escape',char(base)) '\s*\('];
        if ~isempty(regexp(src,pat,'once'))
            blockers(end+1,1) = string(f) + " -> " + legacy(j);
        end
        quoted = ['["' char(39) ']' regexptranslate('escape',char(base)) '\.m["' char(39) ']'];
        if ~isempty(regexp(src,quoted,'once'))
            blockers(end+1,1) = string(f) + " -> " + legacy(j);
        end
    end
end
blockers = unique(blockers);
unused = isempty(blockers);
if unused
    fprintf("Legacy script dependencies: NONE (static scan)\n");
else
    fprintf("Legacy script dependencies found: %d\n",numel(blockers));
    disp(blockers);
end
fprintf("Legacy removal candidates: %d\n",numel(legacy));
end
