# Simulation of a Communication System

## 📌 Overview
This project simulates a communication system that transmits and processes two audio signals. It includes:
- Analysis of audio signals in the time and frequency domains.
- Noise removal using filtering techniques.
- Multiplexing strategies to share a channel with minimal interference.
- Signal recovery using demodulation and filtering.

## 📚 Project Structure
```
📁 Simulation-of-Communication-System
│── 📄 README.md  # Project documentation
│── 📁 src/        # Source code for simulation
│── 📁 audio/      # Audio files used for analysis
│── 📁 filters/    # Filtering algorithms and scripts
│── 📁 results/    # Output and analysis results
│── 📄 report.pdf  # Detailed technical report
```

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

## 🔧 Installation & Usage
### **1. Clone the Repository**
```bash
git clone https://github.com/your-username/Simulation-of-Communication-System.git
cd Simulation-of-Communication-System
```

### **2. Install Dependencies**
```bash
pip install -r requirements.txt
```

### **3. Run the Simulation**
```bash
python src/simulation.py
```

### **4. View Results**
- Processed audio files are saved in the `results/` directory.
- Graphical analysis is available in the `results/plots/` folder.

## 📊 Performance Metrics
The efficiency of filtering and demodulation is measured using the **Signal-to-Interference Ratio (SIR)**. A higher SIR indicates better signal quality. 

| Filter Type                         | SIR (dB) | Best Cutoff Frequency (Hz) |
|--------------------------------------|----------|---------------------------|
| Bessel Low-Pass Filter (1st song)    | 22       | 3200                      |
| Bessel Low-Pass Filter (2nd song)    | 7.3      | 4300                      |
| Narrowband Rejection (1st song)      | 38       | 5537.9                    |
| Narrowband Rejection (2nd song)      | 23       | 5537.9                    |

## 📎 Resources
- [Project Report](report.pdf) (detailed explanation and analysis)
- [Audio & Channel Files](https://drive.google.com/drive/folders/1rXr1-bNQwBm6VB0d-Prug1Dr7H73kD1Z?usp=drive_link)

## 💜 License
This project is licensed under the MIT License.

## 💡 Author
👤 **Dustmamat Bozorkulov**  
🎓 B.Sc. in Electronic and Communications Engineering  
🏢 Politecnico di Torino  
📧 Email: your-email@example.com
