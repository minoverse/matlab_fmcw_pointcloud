function all_pass = run_phase2_regression_suite()
%% Phase 2 실험 전체 회귀 검증
root = fileparts(mfilename("fullpath"));
cd(root);

cases = ["base","2","3a","3b","3c","3c2","3d","4a","4b","4b2","4c","5a","5b","5c","5d","5e"];
passed = false(size(cases));

for k = 1:numel(cases)
    fprintf("\n========== Group %s ==========\n", cases(k));
    try
        rng(1);
        run_phase2(cases(k));
        passed(k) = true;
        fprintf("Group %s: PASS\n", cases(k));
    catch ME
        fprintf("Group %s: FAIL\n", cases(k));
        fprintf("%s\n", getReport(ME,"extended","hyperlinks","off"));
    end
    close all;
end

fprintf("\n========== VALIDATION SUMMARY ==========\n");
for k = 1:numel(cases)
    if passed(k)
        fprintf("Group %s: PASS\n", cases(k));
    else
        fprintf("Group %s: FAIL\n", cases(k));
    end
end
fprintf("Passed: %d / %d\n", nnz(passed), numel(cases));
all_pass = all(passed);
end
