clear all
close all
clc

%% Simulation Parameters
SpS = 64;                  
M = 8;                    
k = log2(M);              
numBits = 2^18 * k;       
numSymbols = numBits / k; 

% Let's define the constellation points (I, Q) explicitly:
constellation_points = [
    -3 -1;  % 000
    -1 -1;  % 010
     1 -1;  % 110
     3 -1;  % 100
    -3  1;  % 001
    -1  1;  % 011
     1  1;  % 111
     3  1   % 101
];

%% Generation of signal and Complex Envelope Signal
bits = randi([0, 1], 1, numBits);
bits_reshaped = reshape(bits, k, numSymbols).'; % Reshape bits into k-bit symbols

% Convert each k-bit sequence to a decimal index
decimal_indices = bi2de(bits_reshaped, 'left-msb') + 1; 

% Look up I and Q values from the constellation
tx_I = constellation_points(decimal_indices, 1);
tx_Q = constellation_points(decimal_indices, 2);

% Combine I and Q into a complex baseband signal
tx_baseband_symbols = tx_I + 1i * tx_Q;

signal_I = repelem(tx_I, SpS); % upsampling
signal_Q = repelem(tx_Q, SpS);
tx_complex_envelope = signal_I + 1i * signal_Q; 

Fs = SpS;
NFFT = 2^nextpow2(length(tx_complex_envelope)); 
Y = fftshift(fft(tx_complex_envelope, NFFT));
deltaFrequency = Fs/(NFFT);
f = (-Fs/2:deltaFrequency:Fs/2-deltaFrequency);

Ps = mean(abs(tx_complex_envelope).^2);
H = fftshift(sinc(f));
%%  Plot the noiseless constellation, both considering only the optimal sampling and the full trajectory
rx_complex_envelope_noiseless =ifft(fft(tx_complex_envelope).*H');
rx_sampled_symbols_noiseless = rx_complex_envelope_noiseless(SpS/2:SpS:end);
figure;
plot(constellation_points(:,1), constellation_points(:,2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
plot(real(rx_sampled_symbols_noiseless), imag(rx_sampled_symbols_noiseless), 'r-', 'LineWidth', 0.5); 
plot(real(rx_sampled_symbols_noiseless), imag(rx_sampled_symbols_noiseless), 'x', 'MarkerSize', 6,'Color', 'g');
grid on;
axis equal;
title('Noiseless Constellation: Optimal Sampled and Full Trajectory');
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
legend('Constellation Points', 'Full Trajectory', 'Optimal Sampled Points');
max_val_noiseless = max(abs(constellation_points(:))) * 1.2; % Dynamic limits
xlim([-max_val_noiseless max_val_noiseless]);
ylim([-max_val_noiseless max_val_noiseless]);
