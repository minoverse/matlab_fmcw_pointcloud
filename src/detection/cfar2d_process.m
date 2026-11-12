
function [cfar_map, noise_map, det_range_idx, det_doppler_idx, stats] = cfar2d_process(RD_power, cfg)

%% CFAR 파라미터
Tr = cfg.Tr;
Td = cfg.Td;

Gr = cfg.Gr;
Gd = cfg.Gd;

Pfa = cfg.Pfa;


%% Training Cell 개수 계산
total_cells = (2*(Tr+Gr)+1) * (2*(Td+Gd)+1);

guard_cells = (2*Gr+1) * (2*Gd+1);

Ntraining = total_cells - guard_cells;


%% 검출 임계값 계수
alpha = Ntraining * (Pfa^(-1/Ntraining) - 1);


%% 검출 및 잡음 맵 초기화
cfar_map = false(size(RD_power));

noise_map = zeros(size(RD_power));

[Nrange, Ndoppler] = size(RD_power);


%% 2D CA-CFAR 검출
for r = Tr+Gr+1 : Nrange-(Tr+Gr)

    for d = Td+Gd+1 : Ndoppler-(Td+Gd)

        % CUT 주변의 윈도우 추출
        window = RD_power(r-(Tr+Gr):r+(Tr+Gr), d-(Td+Gd):d+(Td+Gd));

        training_mask = true(size(window));

        center_r = Tr+Gr+1;
        center_d = Td+Gd+1;

        % Guard Cell과 CUT 제외
        training_mask(center_r-Gr:center_r+Gr, center_d-Gd:center_d+Gd) = false;

        % Training Cell 추출
        training_cells = window(training_mask);

        % 평균 잡음 전력 계산
        noise_power_cfar = mean(training_cells);

        noise_map(r, d) = noise_power_cfar;

        % 검출 임계값 계산
        threshold = alpha * noise_power_cfar;

        CUT = RD_power(r, d);

        % 임계값보다 크면 target으로 검출
        if CUT > threshold
            cfar_map(r, d) = true;
        end

    end

end


%% 검출된 target 인덱스 추출
[det_range_idx, det_doppler_idx] = find(cfar_map);


%% 저장
stats.total_cells = total_cells;
stats.guard_cells = guard_cells;
stats.Ntraining = Ntraining;
stats.alpha = alpha;
stats.Pfa = Pfa;
stats.detected_cells = length(det_range_idx);

end
