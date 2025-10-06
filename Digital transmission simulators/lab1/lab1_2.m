clear all
close all
clc

% Load audio data
[data, Fs] = audioread('Time_10.m4a'); 
N= length(data);
x = data(1:N, 1)'; 
xmax = max(x);
xmin = min(x);
x_mean = mean(x); 
x = x-x_mean; % shift x so that mean value is zero
x_uniform = (xmax - xmin)*rand(1, N) + xmin;

% uniform quantization 
n = 5; 
M = 2^n; % number of quantization levels
delta = (xmax - xmin) / M;
partition = xmin + delta * (1:M-1); 
codebook = xmin + delta * (0.5 : 1 : M - 0.5); % midpoints of the quantization intervals
[indexes, xq] = quantiz(x, partition, codebook);

% tranmission over binary symmetric channel with P_e
Pe_values = logspace(-2, -10, 100);
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
SNR_b = 1./(4.*Pe_values.*((M^2 - 1) / M^2)).*var(x)/var(x_uniform); % SNR for binary sym channel with 
SNR_q = M^2*var(x)/var(x_uniform); % SNR for quantization
SNR_theoretical = 10 * log10(SNR_b.*SNR_q./(SNR_b + SNR_q)); % overall SNR in dB
a = 10 * log10(var(x_uniform)/var(x)) % performance loss

figure(1);
plot(x(1:10000), 'b');
hold on;
plot(xq(1:10000), 'r--');
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
