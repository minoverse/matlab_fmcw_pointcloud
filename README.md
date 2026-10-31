# FMCW Radar DSP with Detection-Centric Sparse Processing

핵심 목표:

- Synthetic ADC → Point Cloud 전체 파이프라인 구현
- Range / Doppler / CFAR / AoA 직접 구현 및 검증
- Dense processing → Detection-Centric Sparse Processing
- 연산량, 메모리, 데이터 이동량 감소
- 절감률과 검출 성능 손실의 trade-off 분석

---

## 1. 전체 구조

```text
Multiple Targets + Noise
        ↓
ADC Cube
        ↓
Range Window
        ↓
Range FFT
        ↓
Static / DC Handling
        ↓
Doppler Window
        ↓
Doppler FFT
        ↓
Range-Doppler Map
        ↓
2D CA-CFAR
        ↓
Peak Grouping
        ↓
Detection List
        ↓
1TX × 4RX ULA
        ↓
AoA / Angle FFT
        ↓
Range + Velocity + Angle + SNR
        ↓
x, y, vr, SNR
        ↓
Point Cloud
        ↓
CSV / Text
```

현재 단계:

Multiple Targets + Noise        ✅
Range FFT + Hann Windowing      ✅
Doppler FFT                     ✅
2D CA-CFAR + Peak Grouping      ✅
Static / DC Handling            ✅
1TX × 4RX AoA                   ✅
Point Cloud + CSV               ✅
Ground-Truth Validation         ✅
Basic Robustness Tests          ✅
Near-Zero Velocity Trade-Off    ✅
Optimization v2                 ← NEXT

---

## 2. ADC Cube

```text
adc_cube(sample, rx, chirp)
```

```text
dim 1 → Fast Time
dim 2 → RX Channel
dim 3 → Slow Time / Chirps
```

처리 방향:

```text
Fast Time   → Range FFT    → Range
Slow Time   → Doppler FFT  → Velocity
RX Channels → Angle FFT    → Angle
```

---

## 3. Baseline Signal Model

Synthetic data:

- Multiple targets
- Range
- Radial velocity
- Amplitude
- Angle
- Complex Gaussian noise

추가 검증 조건:

- Weak targets
- Close-range targets
- Close-velocity targets
- Close-angle targets
- Static clutter
- Near-zero velocity targets
- RX gain / phase mismatch
- SNR sweep

---

## 4. Range / Doppler Processing

### Range

```text
ADC
 ↓
Hann Window
 ↓
Range FFT
 ↓
Range Spectrum
```

검증:

- Range resolution
- Ground Truth Range
- Estimated Range
- Range Error

Range FFT는 Baseline과 Optimization v2 모두 **dense processing**으로 유지.

### Doppler

```text
Range FFT Output
        ↓
Doppler Window
        ↓
Doppler FFT
        ↓
fftshift
        ↓
Range-Doppler Map
```

검증:

- Velocity resolution
- Maximum unambiguous velocity
- Ground Truth Velocity
- Estimated Velocity
- Velocity Error

---

## 5. Static / DC Handling

Strong static clutter가 ROI selection을 방해할 수 있으므로 Slow-Time mean subtraction을 비교 실험.

```text
Range Data
    ↓
Slow-Time Mean
    ↓
Mean Subtraction
```

주의:

```text
v ≈ 0 target
Stopped vehicle
Guardrail
Tangential motion
```

도 함께 억제될 수 있음.

따라서:

```text
With Clutter Removal
vs
Without Clutter Removal
```

을 비교.

Static removal은 항상 적용하는 전처리가 아니라 **trade-off 항목**으로 처리.

---

## 6. 2D CA-CFAR + Peak Grouping

```text
Training Cells
      +
Guard Cells
      +
CUT
      ↓
Noise Estimate
      ↓
Threshold
      ↓
Detection
```

CFAR output:

```text
Range Bin
Doppler Bin
Detection Power
```

한 target이 여러 cell로 검출될 수 있으므로:

```text
CFAR Detections
        ↓
Local Peak
        ↓
Peak Grouping
        ↓
Detection List
```

---

## 7. 1TX × 4RX AoA

Baseline:

```text
1 TX
×
4 RX ULA
```

Target angle에 따라 RX antenna 사이에 spatial phase progression 발생.

```text
RX1
RX2
RX3
RX4
 ↓
Spatial Phase Difference
 ↓
Angle FFT
 ↓
Azimuth Estimate
```

검증:

- Ground Truth Angle
- Estimated Angle
- Angle Error
- Close-angle target separation

---

## 8. Point Cloud

각 detection에서:

```text
Range
Velocity
Angle
SNR
```

추정.

좌표계:

```text
x → Forward
y → Lateral
```

```text
Range + Angle
      ↓
Cartesian Conversion
      ↓
x, y
```

최종 Point Cloud:

```text
x
y
range
radial_velocity
angle
SNR
```

---

## 9. Baseline Validation

Optimization 전에 Baseline 정확성을 먼저 확인.

### 이론값

- Range resolution
- Velocity resolution
- Maximum unambiguous velocity
- FFT bin spacing

### Ground Truth

- Range Error
- Velocity Error
- Angle Error

### Robustness

- Weak targets
- Multiple targets
- Close targets
- Static clutter
- Near-zero velocity
- RX mismatch
- SNR sweep

---

# Optimization v2

## 10. Embedded Sparse DSP Optimization
```text
Dense Range Processing
        ↓
Lightweight Early Gating
        ↓
Candidate Range Bins
        ↓
Selective Downstream Processing
        ↓
Sparse Point Cloud
```
```text
목표:

Unnecessary Processing ↓
Compute ↓
Memory Footprint ↓
Data Movement ↓
```
동시에:

Detection Performance 유지
Weak-Target Miss 최소화
Range / Velocity / Angle Accuracy 유지

## 11. Optimization Architecture

```text
Range FFT                    ← Dense
    ↓
Static / DC Handling
    ↓
Adaptive ROI Selection
    ↓
Candidate Range Bins
    ↓
Selective Doppler FFT
    ↓
ROI-Only CFAR
    ↓
Peak Grouping
    ↓
Detection List
    ↓
Detection-Based AoA
    ↓
Sparse Point Cloud
```

주요 절감 대상:

- Doppler FFT
- CFAR processing
- Buffer usage
- Data movement
- Output data size

---

## 12. ROI Candidate Selection

Range FFT 이후 low-cost candidate 생성.

```text
Range FFT
    ↓
Static / DC Handling
    ↓
Noncoherent Energy
    ↓
Noise Floor Estimate
    ↓
Median × Scale
    ↓
Candidate Range Bins
```

ROI threshold는 conservative하게 설정.

이유:

```text
Pre-Doppler ROI
→ Doppler integration gain 없음
→ Weak target miss 가능
```

핵심 평가:

```text
Processing Reduction
        vs
Missed Target Rate
```

절감률이 작더라도 weak target을 유지하는 것이 우선.

---

## 13. Selective Doppler + ROI CFAR

### Baseline

```text
All Range Bins
      ↓
Doppler FFT
      ↓
Full CFAR Search
```

### Sparse

```text
Candidate Range Bins
        ↓
Selective Doppler FFT
        ↓
ROI-Only CFAR
```

CFAR ROI는 edge detection을 위해:

```text
Candidate ROI
      ↓
Training + Guard Margin
      ↓
ROI-Only CFAR
```

으로 확장.

Margin cell도 실제 processing cost에 포함.

---

## 14. Fair Comparison

### Dense Reference

```text
Full Range
→ Full Doppler
→ Dense Downstream
```

### Practical Baseline

```text
Full Range FFT
→ Full Doppler FFT
→ Full CFAR
→ Detection List
→ AoA
→ Point Cloud
```

### Sparse Optimization

```text
Full Range FFT
→ ROI Selection
→ Selective Doppler FFT
→ ROI-Only CFAR
→ Detection List
→ AoA
→ Sparse Point Cloud
```

실제 핵심 비교:

```text
Practical Baseline
        vs
Sparse Optimization
```

Detection-based AoA 자체는 Optimization 성과로 과장하지 않음.

---

## 15. 평가 지표

### Radar Performance

```text
Pd
Pfa
Missed Target Rate
Range RMSE
Velocity RMSE
Angle RMSE
```

### Processing Cost

```text
Selected Range Bins
Doppler FFT Count
CFAR CUT Count
Operation Count
Buffer Size
Estimated Bytes Moved
Output Data Size
```

Pfa는 단일 frame이 아니라 noise-only Monte Carlo로 측정.

```text
Noise-Only Frames
      ↓
Repeated CFAR
      ↓
False Detections / Tested CUTs
      ↓
Measured Pfa
```

---

## 16. Worst-Case / Reliability

Sparse processing은 scene이 sparse할수록 유리.

```text
More Targets / Clutter
        ↓
More ROI Bins
        ↓
More Doppler FFT
        ↓
More CFAR Processing
        ↓
Higher Workload
```

따라서:

```text
Typical Case
vs
Worst Case
```

를 함께 측정.

핵심 관점:

> 평균 연산량뿐 아니라 scene이 복잡해졌을 때 workload upper bound도 확인.

---

## 17. 핵심 Trade-offs

### Static Removal

```text
Clutter 제거 ↑
Near-zero target sensitivity ↓
```

### ROI Threshold

```text
Low Threshold
→ Miss ↓
→ Saving ↓

High Threshold
→ Saving ↑
→ Miss ↑
```

### ROI-CFAR

```text
Smaller ROI
→ CFAR Cost ↓

Training + Guard Margin
→ 실제 Saving 감소
```

### Hardware

Algorithmic operation 감소가 실제 power / latency 감소와 동일하지는 않음.

현재 평가는:

```text
Operations
Processed Cells
Memory
Estimated Data Movement
```

중심.

---

## 18. Optional — Fixed-Point

```text
Floating Point
      ↓
8 / 12 / 16-bit
      ↓
Quantization
      ↓
FFT Scaling
      ↓
Overflow / Saturation
```

비교:

- Memory
- Quantization noise
- SNR loss
- Pd
- RMSE
- Overflow

---

## 19. Advanced Extension — 3TX × 4RX TDM-MIMO

Sparse Optimization과 별도 확장.

```text
1TX × 4RX
      ↓
3TX × 4RX TDM-MIMO
      ↓
TX Sequencing
      ↓
Virtual Antenna Geometry
      ↓
Doppler Compensation
      ↓
AoA
      ↓
Gain / Phase Mismatch
      ↓
Calibration
```

주요 관심:

- Virtual array
- TDM Doppler compensation
- Channel mismatch
- Calibration
- AoA workload scaling

---

## 20. Main Contribution

```text
Radar DSP
+
Compute-on-Demand
+
Event-Driven Thinking
+
Memory Awareness
+
Data Movement Awareness
+
Worst-Case Reliability
```

핵심 질문:

> **유용한 Radar 정보를 잃기 시작하기 전까지 processing을 얼마나 줄일 수 있는가?**

Synthetic data에서 먼저 검증하고, 이후 real radar data에서 동일한 경계를 다시 확인.

---

## References

### Texas Instruments

- TI mmWave SDK  
  https://www.ti.com/tool/MMWAVE-SDK

- TI mmWave MCU+ SDK  
  https://www.ti.com/tool/MMWAVE-MCUPLUS-SDK

- TI mmWave-L-SDK  
  https://www.ti.com/tool/MMWAVE-L-SDK

- TI AWR1843  
  https://www.ti.com/product/AWR1843

- TI AWR1843 Datasheet  
  https://www.ti.com/lit/ds/symlink/awr1843.pdf

### MathWorks

- Range-Doppler Response  
  https://www.mathworks.com/help/phased/ug/range-doppler-response.html

- Range and Doppler Estimation  
  https://www.mathworks.com/help/phased/range-doppler-estimation.html

- 2-D CFAR Detector  
  https://www.mathworks.com/help/phased/ref/phased.cfardetector2d-system-object.html

---


→ Optimization v2
```
