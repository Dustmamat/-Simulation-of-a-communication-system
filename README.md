# Simulation of a Communication System

## 📌 Overview
This project simulates a communication system that transmits and processes two audio signals. It includes:
- Analysis of audio signals in the time and frequency domains.
- Noise removal using filtering techniques.
- Multiplexing strategies to share a channel with minimal interference.
- Signal recovery using demodulation and filtering.


## 🎯 Key Features
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


## 📊 Performance Metrics
The efficiency of filtering and demodulation is measured using the **Signal-to-Interference Ratio (SIR)**. A higher SIR indicates better signal quality. 

| Filter Type                         | SIR (dB) | Best Cutoff Frequency (Hz) |
|--------------------------------------|----------|---------------------------|
| Bessel Low-Pass Filter (1st song)    | 22       | 3200                      |
| Bessel Low-Pass Filter (2nd song)    | 7.3      | 4300                      |
| Narrowband Rejection (1st song)      | 38       | 5537.9                    |
| Narrowband Rejection (2nd song)      | 23       | 5537.9                    |

## 📎 Resources
- [Audio & Channel Files](https://drive.google.com/drive/folders/1rXr1-bNQwBm6VB0d-Prug1Dr7H73kD1Z?usp=drive_link)

## 💜 License
This project is licensed under the MIT License.

## 💡 Author
👤 **Dustmamat Bozorkulov**  
🎓 B.Sc. in Electronic and Communications Engineering  
🏢 Politecnico di Torino  
📧 Email: s307576@studenti.polito.it
