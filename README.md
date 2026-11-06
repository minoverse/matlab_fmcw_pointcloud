
# FMCW Radar DSP — Detection-Centric Sparse Processing

## 1. 프로젝트 개요

MATLAB 기반 FMCW 레이더 신호처리 파이프라인 구현 및 **표적 검출 영역을 선택적으로 처리하는 Sparse Processing 최적화** 프로젝트.

**목표:** 불필요한 연산 감소와 함께 검출 성능, 메모리 접근 비용, 처리량 제한 및 표적 누락 위험 검증.

### 시스템 구성

- 레이더: 77 GHz FMCW, 대역폭 1 GHz, 64 Chirps, 1,000 Samples/Chirp
- 안테나: 4 RX (2×2 Planar Array), 방위각·고도각 추정
- 출력: 거리·속도 검출 및 3D Point Cloud

**처리 흐름**

ADC → Range FFT → ROI 선택 → Selective Doppler FFT → CFAR → 각도 추정 → 3D Point Cloud

## 2. 주요 최적화 결과

| 검증 항목 | 결과 |
|---|---|
| Doppler FFT 횟수 | 2,000 → 176 (**91.2% 감소**) |
| 전체 FFT 연산량 | **21.10% 감소** (구조적 상한 23.14%) |
| CFAR 검사 셀 | 25,920 → 216 (**99.17% 감소**) |
| 모델 기반 데이터 접근량 | **29.46% 감소** |
| Dense / Sparse 수치 일치 | 공통 검출 표적의 추정 결과 일치 (**PASS**) |

*단일 검증 장면 기준. 전체 FFT 연산량은 ROI 관리 오버헤드 제외, 데이터 접근량은 Dense와 Sparse에 동일한 Integral-Image CFAR 모델을 적용한 추정치.*

### 검출 성능 및 한계

- **SNR 20 dB, ROI 3 dB:** Dense / Sparse 모두 60/60 검출
- **SNR 20 dB, ROI 6 dB:** Dense 60/60, Sparse 30/60 — 30건 추가 누락
- **저 SNR:** 약한 표적이 ROI 선택 단계에서 제외되는 문제 확인
- **밀집 장면:** FFT 예산 준수, Dense 검출 48건 중 6건만 보존

**핵심 발견:** 연산량 절감과 검출 성능 사이의 Trade-off 존재. ROI 임계값 및 처리 예산을 검출 손실과 함께 고려할 필요.

Dense / Sparse 공통 검출의 수치 일치는 회귀 검증 결과이며, 정확도 향상이나 모든 조건에서의 표적 보존을 의미하지 않음.

## 3. 검증 결과 시각화

### 3.1 3D Point Cloud — 전체 DSP 구현 결과

거리·속도·각도 추정부터 3차원 표적 위치 생성까지의 처리 결과.

![3D Point Cloud](plots/06_pointcloud/pointcloud.png)

### 3.2 ROI 임계값 — 검출 성능과 연산량의 Trade-off

ROI 선택 조건에 따른 검출 손실과 Doppler 처리량의 관계.

![ROI Threshold vs Detection and Workload](output/phase2/04_roi_threshold_sweep_miss_rate/roi_threshold_miss_rate_vs_doppler_workload.png)

### 3.3 약한 표적 검출 — CFAR 보정 효과와 한계

CFAR 임계값 보정에 따른 약한 표적 검출 특성 비교.

![Calibrated CFAR Weak Target Trade-off](output/phase2/07_calibrated_cfar_detection_validation/calibrated_cfar_weak_target_tradeoff.png)

## 4. 프로젝트 구성 및 실행

```text
config.m                     # 레이더 설정 및 Ground Truth
main.m                       # 기본 DSP 및 3D Point Cloud

src/
  processing/                # Range / Doppler FFT
  detection/                 # CFAR / Peak Grouping
  angle/                     # 방위각 / 고도각 추정
  output/                    # Point Cloud / 결과 검증
  optimization/              # ROI / Top-K / Sparse Processing

phase2_group*.m              # 단계별 최적화 및 검증 실험
phase2_group5d_main.m        # CFAR 데이터 접근 비용 비교
phase2_group5e_main.m        # Dense / Sparse 정확도 검증

tests/                       # 회귀 테스트
output/phase2/               # 실험 결과 (CSV / TXT / PNG)
plots/                       # DSP 및 Point Cloud 시각화
```

**MATLAB 실행 — 프로젝트 루트(`config.m` 위치)**

```matlab
addpath(genpath("src"));

main
phase2_group5d_main
phase2_group5e_main
```

**주요 검증 자료**

- [Group 5D — CFAR 데이터 접근 비용 검증](output/phase2/15_fair_cfar_cost_validation/fair_cfar_cost_validation.txt)
- [Group 5E — Dense / Sparse 정확도 검증](output/phase2/16_dense_sparse_accuracy_validation/dense_sparse_accuracy_validation.txt)
- [Group 5C — Monte Carlo 통계 검증](output/phase2/14_final_statistical_validation/final_statistical_validation.txt)
- [자원 사용량 비교](output/phase2/12_unified_resource_counters/key_resource_comparison.csv)

## 5. Engineering Perspective

**Local Optimization ≠ System-Level Optimization**

FFT 연산량 감소뿐 아니라 다음 항목을 함께 고려한 검증.

- CFAR 처리를 위한 ROI 확장 비용
- 전체 연산량에서 Dense Range FFT가 차지하는 비중
- 저 SNR 및 밀집 장면에서의 표적 누락
- 메모리 접근 비용 및 제한된 처리 예산

**핵심 원칙:** 검출 성능과 자원 제약을 고려해 필요한 연산만 선택적으로 수행.

**검증 범위:** MATLAB 시뮬레이션 및 알고리즘 검증 완료. 실제 Radar ADC 데이터, 하드웨어 실행시간, 물리적 메모리 사용량, 버스 트래픽 및 전력 소비는 미측정.
