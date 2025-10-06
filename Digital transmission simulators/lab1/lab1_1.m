clear all
close all
clc

N = 100000; % number of samples
xmax = 1;
xmin = -1;
x = (xmax - xmin)*rand(1, N) + xmin; % signal  with uniform pdf in range [-1, 1]

% uniform quantization 
n = 5; % number of bits
M = 2^n; % number of quantization levels
delta = (xmax - xmin) / M;
partition = xmin + delta * (1:M-1); 
codebook = xmin + delta * (0.5 : 1 : M - 0.5); % midpoints of the quantization intervals
Pe_values = logspace(-2, -10, 100);
[indexes, xq] = quantiz(x, partition, codebook);


% tranmission over binary symmetric channel with P_e
bits_tx= de2bi(indexes,n);
SNR_simulation = zeros(size(Pe_values));

for k = 1:length(Pe_values)
    pe = Pe_values(k);
    bits_rx = bsc(bits_tx, pe);
    indexes_out= bi2de(bits_rx);
    x_out = codebook(indexes_out + 1);
    noise = x - x_out;
    SNR_simulation(k) = 10 * log10(var(x) / var(noise)); % in 
end

% theoretical estimation of SNR 
SNR_b = 1./(4.*Pe_values.*((M^2 - 1) / M^2));
SNR_q = M^2;
SNR_theoretical = 10 * log10(SNR_b.*SNR_q./(SNR_b + SNR_q)); % in dB


figure(1);
plot(x(1:100), 'b');
hold on;
plot(xq(1:100), '--');
legend('Original Signal', 'Quantized Signal');
title('Original vs Quantized Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

figure(2);
histogram(x, 100, 'Normalization', 'pdf');
xlabel('Amplitude');
ylabel('Probability Density');
grid on;

figure(3);
semilogx(Pe_values, SNR_simulation, 'b-', 'LineWidth', 2);
hold on;
semilogx(Pe_values, SNR_theoretical, 'r--', 'LineWidth', 2);
xlabel('Bit Error Probability (P_e)');
ylabel('SNR (dB)');
legend('Simulation', 'Theoretical');
grid on;
