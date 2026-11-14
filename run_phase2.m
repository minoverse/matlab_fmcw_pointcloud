function run_phase2(group)

%% 프로젝트 루트 

root = fileparts(mfilename("fullpath"));
cd(root);


%% 모듈 및 test 경로 등록

addpath(root);
addpath(genpath(fullfile(root, "src")));
addpath(fullfile(root, "tests"));

exp_dir = fullfile(root, "experiments", "phase2");
addpath(exp_dir);


%% test 선택

if string(group) == "base"

    phase2_main;

else

    script_name = "phase2_group" + string(group) + "_main";

    script_file = fullfile(exp_dir, script_name + ".m");

    assert(isfile(script_file), "Experiment not found");

    assert(strcmp(which(char(script_name)), char(script_file)));

    feval(char(script_name));

end

end
