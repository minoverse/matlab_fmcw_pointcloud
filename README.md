
# FMCW Radar DSP with Detection-Centric Sparse Processing

MATLAB 기반 FMCW 레이더 DSP 구현 및 **검출 손실과 Worst-Case 처리량을 고려한 Sparse Processing 최적화** 프로젝트.

**목표:** 후단 연산량을 줄이면서 표적 검출 성능, 숨겨진 처리 비용, 고정된 FFT 예산에서의 실패 조건을 함께 검증.

## 1. 시스템

- **레이더:** 77 GHz FMCW, Bandwidth 1 GHz, 1,000 Samples/Chirp, 64 Chirps
- **안테나:** 4 RX (2×2 Planar Array), Azimuth/Elevation 추정
- **출력:** Range–Doppler Detection, Peak Grouping, 3D Point Cloud
- **처리:** `ADC → Dense Range FFT → ROI Energy Gating → Selective Doppler FFT → ROI-Only CFAR → Peak Grouping → AoA → 3D Point Cloud`
- Static/DC 제거는 Near-Zero Velocity 표적 손실 가능성을 고려해 선택적으로 적용

## 2. 핵심 결과

| 항목 | 결과 및 제한 조건 |
|---|---|
| Doppler FFT | Fair Dense **2,000 → 176** (Group 5A 사례, **91.20% 감소**); Group 5C 평균 **176.98** |
| Range + Doppler FFT `N log₂N` | **21.10% 감소**; Doppler만 제거할 때 구조적 상한 **23.14%** (Range FFT 비중 76.86%) |
| CFAR CUT | **25,920 → 216** (**99.17% 감소**, Group 5A 단일 사례); Group 5E 20 dB / ROI 3 dB 평균 **232.2 CUT** (약 **99.10% 감소**) |
| 논리적 버퍼 합계 *(모델)* | **10.284 → 8.008 MiB** (**22.13% 감소**) |
| 데이터 접근량 *(공정 Integral-Image CFAR 모델)* | **14.475 → 10.210 MiB** (**29.46% 감소**, Group 5D **단일 장면**, PASS) |
| Naive CFAR 접근 모델 *(참고)* | 54.434 → 10.465 MiB (80.77%); 구현 방식에 민감하여 대표 절감률로 사용하지 않음 |
| Group 5C 약한 표적, SNR 20 dB, ROI 3 dB | Sparse **200/200** local hits (amplitude 0.005, CFAR α=4.8) |
| 동일 조건, ROI 4 dB / 5–6 dB | Sparse **59/200 / 0/200** local hits |
| SNR 10 dB, ROI 3 dB *(Group 5C)* | Dense **182/200**, Sparse **0/200** local hits |
| 16 Targets 밀집 장면, FFT Cap 200 | Cap 준수 **PASS**, Dense 검출 48건 중 **6건만 보존** (3 Trials 합계) |
| Group 5E 실제 Sparse FFT/CFAR/AoA | 20 dB ROI 3 dB에서 Dense/Sparse **각 60/60**, 공통 표적 추정값 최대 차이 **0** |
| CFAR Pfa Calibration | 보정 α=4.8, holdout **78 false alarms / 7,776,000 CUTs** |

**Group 5E 절대 RMSE — 20 dB, ROI 3 dB, 공통 검출 60건 (30 Trials, 원래 CFAR α≈11.825)**

| 경로 | Range [m] | Velocity [m/s] | Azimuth [°] | Elevation [°] | 누락 |
|---|---:|---:|---:|---:|---:|
| Dense | 0.0354 | 0.2053 | 1.2407 | 1.2968 | 0/60 |
| Actual Sparse | 0.0354 | 0.2053 | 1.2407 | 1.2968 | 0/60 |

**5E RMSE의 의미:** Dense와 Sparse의 공통 검출에서 추정값이 동일한 것은 **수치적 동등성(Regression) 확인**이지, 잡음에 대한 정확도 개선의 증거가 아님. 같은 강한 표적의 오차가 여러 SNR에서 일정하여 **FFT Bin 양자화 오차가 지배적**인 결과로 해석되며, Noise Sensitivity 평가는 별도로 필요함. RMSE는 **검출된 표적에 대한 조건부 값**이므로 반드시 누락 수와 함께 해석.

**누락 해석:** Group 5E의 20 dB ROI 6 dB에서 Dense **60/60**, Sparse **30/60**으로 **30건 추가 누락**. 공통 30건의 RMSE 일치가 누락된 30건의 정확성을 뜻하지 않음. **10 dB ROI 3 dB에서는 Dense가 이미 약한 표적을 놓쳐 양쪽 모두 30/60**이므로 Sparse 추가 누락 0건을 **Sparse 안전성의 근거로 사용할 수 없음**. Group 5C는 다른 CFAR 계수·Trial/Seed를 사용한 별도 검증임.

**해석상의 한계:** 21.10%는 FFT 연산량 지표이며 Gating Energy/Median, 정렬, 스케줄링, 버퍼 관리 등 오버헤드를 제외. 29.46%는 **Group 5D 단일 장면의 논리적 접근 모델** 결과로, 실제 메모리 버스 트래픽·전력·지연시간 측정치가 아님. Sparse ROI Energy Read **1.953 MiB**를 Range FFT 출력의 추가 읽기로 가정했으며, 계산을 FFT 출력 생성과 융합하면 이 비용이 달라질 수 있음. FFT Cap PASS가 표적 보존을 보장하지 않음.

**Group 5C 통계 범위:** 1,600개 독립 Noise Frame에서 Dense CFAR 판단에 ROI Mask를 적용한 결과이며 Sparse FFT/CFAR를 1,600회 재실행한 것은 아님. ±1 Range/Doppler Bin의 Local Hit는 저 SNR의 우연한 Noise Hit를 포함할 수 있음. 실제 Sparse 연산의 수치 일치는 Group 4/5D/5E에서 별도 확인.

## 3. 검증 그래프

![FFT workload and structural limit](docs/figures/fft_workload.png)

![Weak-target detection versus ROI threshold](docs/figures/weak_target_threshold.png)

## 4. 프로젝트 구성 및 실행

```text
config.m                       # FMCW 설정 / Ground Truth
main.m                         # Dense Baseline + 3D Point Cloud
src/processing/                # Range/Doppler FFT
src/detection/                 # CFAR / Peak Grouping
src/angle/                     # 2×2 Planar Array AoA
src/output/                    # Point Cloud / Ground Truth 검증
src/optimization/              # ROI / Top-K / Block / Sparse CFAR
phase2_group*.m                # 단계별 최적화 및 검증 실험
  phase2_group5d_main.m        # 공정 CFAR 접근량 모델 검증
  phase2_group5e_main.m        # 실제 Sparse FFT/CFAR/AoA / RMSE 검증
tests/                         # Regression 테스트
output/phase2/06_*/            # CFAR Pfa Calibration
output/phase2/08_*/ ~ 14_*/   # Block, Budget, Resource, Monte Carlo
output/phase2/15_*/ ~ 16_*/   # Fair CFAR Cost / Dense-Sparse Accuracy
docs/                          # 검증 보고서 및 그래프
```

**MATLAB: 프로젝트 루트(`config.m` 위치)에서 실행**

```matlab
addpath(genpath("src"));
phase2_group5d_main
phase2_group5e_main
```

Group 5C는 장시간 수행될 수 있으며 동일 설정의 체크포인트가 있으면 완료된 장면을 재사용함.

## 5. Engineering Perspective

**Local Optimization ≠ System-Level Optimization.** 연산량 감소만으로 최적화 성공을 판단하지 않고 CFAR Margin, Dense Range FFT, 저 SNR 누락, 메모리 접근 모델의 공정성, 밀집 장면의 예산 상한을 함께 평가.

**핵심 원칙:** 검출 성능과 자원 제약을 만족하는 범위에서 필요한 계산만 수행.

**검증 범위:** MATLAB 시뮬레이션 및 알고리즘 검증 완료(Group 5D/5E PASS). 실제 Radar Data, MCU/SoC 실행시간, 물리적 Peak RAM, Bus Traffic 및 Frame당 Energy는 미측정.

