# Digital Communication Systems Simulation

## 📌 Overview
This repository contains a comprehensive series of digital communication system simulations developed as part of the **Digital Transmission**, **Signals and Systems** courses at **Politecnico di Torino**. 

---

##  Projects

### **Project 1: Simulation of a Communication System**

####  Overview
This project simulates a communication system that transmits and processes two audio signals. It includes:
- Analysis of audio signals in the time and frequency domains.
- Noise removal using filtering techniques.
- Multiplexing strategies to share a channel with minimal interference.
- Signal recovery using demodulation and filtering.

####  Key Features
1. **Audio Signal Analysis**:
   - Time and frequency domain visualization.
   - Bandwidth calculation.

2. **Disturbance Removal**:
   - **Bessel Low-Pass Filter**: Removes high-frequency noise while preserving signal integrity.
   - **Narrow Frequency Rejection Filter**: Targets specific noise frequencies without affecting the desired signal.

3. **Channel Sharing via Multiplexing**:
   - Amplitude modulation (AM) used to transmit two audio signals in parallel.
   - Selection of optimal carrier frequencies to prevent interference.

4. **Signal Recovery**:
   - Demodulation using a local carrier and filtering.
   - Performance comparison of different low-pass filters.

####  Performance Metrics
The efficiency of filtering and demodulation is measured using the **Signal-to-Interference Ratio (SIR)**. A higher SIR indicates better signal quality.

| Filter Type                         | SIR (dB) | Best Cutoff Frequency (Hz) |
|--------------------------------------|----------|---------------------------|
| Bessel Low-Pass Filter (1st song)    | 22       | 3200                      |
| Bessel Low-Pass Filter (2nd song)    | 7.3      | 4300                      |
| Narrowband Rejection (1st song)      | 38       | 5537.9                    |
| Narrowband Rejection (2nd song)      | 23       | 5537.9                    |

#### 📎 Resources
- [Project 1 Report](https://drive.google.com/file/d/1YnaELcqHel7hu-QJDJIUZ4qz2SyKdq74/view?usp=drive_link)
- [Audio & Channel Files](https://drive.google.com/drive/folders/1rXr1-bNQwBm6VB0d-Prug1Dr7H73kD1Z?usp=drive_link)

---

### **Project 2: Lab 1 - Analog Signal Transmission over Digital Channel**

####  Overview
This project develops a numerical simulator to evaluate the performance of transmitting analog signals over a digital channel. The simulation analyzes the impact of quantization and Binary Symmetric Channel (BSC) characteristics on signal integrity.

####  Key Features
- **Quantization Analysis**: Fixed 5-bit quantization with performance evaluation
- **BSC Transmission**: Bit error probability ranging from 10⁻² to 10⁻¹⁰
- **Real Signal Testing**: Application to real music signals with non-uniform distributions
- **Huffman Coding**: Data compression technique for improved transmission efficiency

####  Key Results
- **SNR Performance**:
  - Quantization-limited SNR: ~30 dB (5 bits)
  - Channel-limited regime follows SNR ∝ 1/(4Pₑ)
- **Music Signal**: 10.26 dB degradation compared to uniform signal
- **Huffman Coding**: Achieved coding gain of ~1.4
- **Source Entropy**: H = 3.55 bits/symbol

####  Theoretical Foundation
Total SNR for uniform quantization: SNR = M² / [1 + 4Pₑ(M² - 1)]
where M = 2ⁿ (n = number of bits), Pₑ = bit error probability

#### 📎 Resources
- [Lab 1 Report](https://drive.google.com/file/d/1aIeeNznBSbLCcnG8a2UcJuFQMn66d8NY/view?usp=drive_link)

---

### **Project 3: Lab 2 - PAM Modulation with Receiver Filtering Strategies**

####  Overview
This project investigates the performance of digital communication systems using Pulse Amplitude Modulation (PAM-2 and PAM-4) under various receiver filtering strategies. The study compares matched filters with single-pole filters and analyzes system resilience to threshold variations.

####  Key Features
1. **Matched Filter Performance**:
   - Optimal receiver filtering for PAM-2 and PAM-4
   - BER vs. Eb/N0 analysis
   - Eye diagram visualization

2. **Single-Pole Filter Optimization**:
   - Cutoff frequency optimization to minimize Eb/N0 at BER = 10⁻⁴
   - Implementation penalty quantification
   - ISI (Intersymbol Interference) analysis

3. **Threshold Resilience**:
   - Decision threshold variation analysis
   - Maximum allowable threshold error for 0.5 dB penalty

####  Key Results

**Single-Pole Filter Performance:**

| Parameter | PAM-2 | PAM-4 |
|-----------|-------|-------|
| Optimal Cutoff Frequency (normalized) | 0.4 | 0.3 |
| Min Eb/N0 at BER = 10⁻⁴ (dB) | 10.2 | 14.9 |
| Implementation Penalty (dB) | 1.8 | 2.6 |
| Samples per Symbol (SpS) | 8 | 8 |

**Threshold Resilience (PAM-2):**
- Maximum threshold deviation: ±0.08 for 0.5 dB penalty
- Symmetric performance degradation around optimal threshold

####  Simulation Parameters
- **Number of Bits**: 2²⁰ (1,048,576 bits)
- **Samples per Symbol**: 8 (default), tested up to 64
- **Target BER**: 10⁻⁴
- **Modulation**: NRZ rectangular pulses
- **Gray Coding**: Used for PAM-4 to minimize bit errors

#### 📎 Resources
- [Lab 2 Report](https://drive.google.com/file/d/1Jyt5kT7tbsgghHO823bS6_LYcH6Gg-dy/view?usp=drive_link)

---

### **Project 4: Lab 3 - 8-QAM Digital Transmission System**

####  Overview
This project presents a comprehensive analysis of a numerical simulator for an 8-QAM (Quadrature Amplitude Modulation) digital transmission system using NRZ (Non-Return to Zero) pulses. The study compares six different 8-QAM constellation types to evaluate their performance characteristics.

####  Key Features
- **Multiple 8-QAM Constellations**: Rectangular, Circular (8-PSK), Square, SP, X, and Star
- **High Statistical Reliability**: 2¹⁸ (262,144) symbols transmitted
- **Matched Filtering**: 16 samples/symbol for optimal signal detection
- **Nearest Neighbor Decoding**: Euclidean distance-based demodulation
- **Complex Envelope Analysis**: I/Q component representation

####  Constellation Types Analyzed

1. **Rectangular 8-QAM**: Standard grid arrangement
2. **Circular 8-QAM (8-PSK)**: Equal amplitude, phase-based encoding
3. **Square 8-QAM**: Square grid configuration
4. **SP 8-QAM**: Optimized spacing variant
5. **X 8-QAM**: X-shaped constellation
6. **Star 8-QAM**: Star-pattern arrangement

####  Performance Comparison

**Eb/N0 Required at BER = 10⁻⁴:**

| Modulation Type | Eb/N0 (dB) | Performance |
|-----------------|------------|-------------|
| **SP 8-QAM** | **11** | **Best** |
| **Star 8-QAM** | **11** | **Best** |
| Rectangular 8-QAM | 12 | Good |
| Circular 8-QAM | 12 | Good |
| Square 8-QAM | 12 | Good |
| X 8-QAM | 14 | Worst |

**Key Findings:**
- **Best Performers**: SP and Star 8-QAM (3 dB better than X 8-QAM)
- **Performance Factor**: Minimum Euclidean distance between constellation points
- **Optimal Sampling**: Consistently at SpS/2 for all constellation types

####  Theoretical Background
- **Complex Envelope**: sbb(t) = I(t) + jQ(t)
- **Passband Signal**: s(t) = Re{sbb(t) · e^(j2πfct)}
- **Orthogonal Components**: I(t)cos(2πfct) - Q(t)sin(2πfct)

#### 📎 Resources
- [Lab 3 Report](https://drive.google.com/file/d/1rLrFscvd7CEjzRViykXfCbBnf0lzqh_Q/view?usp=drive_link)
