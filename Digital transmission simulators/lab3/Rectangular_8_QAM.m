clear all
close all
clc

%% Simulation Parameters
SpS = 16;                 % Samples per Symbol 
M = 8;                    
k = log2(M);              % Number of bits per symbol (k=3 for 8-QAM)
numBits = 2^18 * k;       % number of bits 
numSymbols = numBits / k; % Total number of symbols

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
% Map bits to symbols (I and Q values)
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

%% 1. Plot the spectrum of the generated signal
figure;
Fs = SpS;
NFFT = 2^nextpow2(length(tx_complex_envelope)); 
Y = fftshift(fft(tx_complex_envelope, NFFT));
deltaFrequency = Fs/(NFFT);
f = (-Fs/2:deltaFrequency:Fs/2-deltaFrequency);

plot(f, 20*log10(abs(Y)));
xlabel('Normalized Frequency (Hz/Symbol_Rate)');
ylabel('Magnitude (dB)');
grid on;

%% 2. Plot the BER vs. Eb/N0 curve
EbN0_dB_range = 2:1:15; 
BER = zeros(1, length(EbN0_dB_range));
theor_BER = zeros(1, length(EbN0_dB_range));
Ps = mean(abs(tx_complex_envelope).^2);
H = fftshift(sinc(f));

for i = 1:length(EbN0_dB_range)
    EbN0_dB = EbN0_dB_range(i);
    EbN0_linear = 10^(EbN0_dB / 10);

    sigma_n_squared = Ps*SpS./(2*k*EbN0_linear);
    sigma_n = sqrt(sigma_n_squared); 

    % Add AWGN to the complex envelope signal
    noise_I = sigma_n * randn(size(tx_complex_envelope));
    noise_Q = sigma_n * randn(size(tx_complex_envelope));
    
    rx_complex_envelope =ifft(fft(tx_complex_envelope + (noise_I + 1i * noise_Q)).*H');
    rx_sampled_symbols = rx_complex_envelope(SpS/2:SpS:end); % to choose sampling instant by using eyediagram 

    % Demodulation: Nearest Neighbor Decoding
    decoded_bits = zeros(1, numBits); 
    
    for s_idx = 1:numSymbols
        received_point = rx_sampled_symbols(s_idx);
        
        % Calculate Euclidean distances from received_point to all constellation points
        distances = abs(received_point - (constellation_points(:,1) + 1i * constellation_points(:,2)));
        
        % Find the index of the closest constellation point
        [~, min_idx] = min(distances);
        
        % Convert the index back to a k-bit sequence
        decoded_k_bits = de2bi(min_idx - 1, k, 'left-msb');
        
        start_idx = (s_idx - 1) * k + 1;
        end_idx = s_idx * k;
        decoded_bits(start_idx:end_idx) = decoded_k_bits;
    end
    
    % Calculate Bit Error Rate (BER)
    num_errors = sum(xor(bits, decoded_bits));
    BER(i) = num_errors / numBits;
    theor_BER(i) =(1/k) * 2 * (1 - 1/sqrt(M)) * erfc(sqrt( (3 * k) / (2 * (M-1)) * EbN0_linear ));
    
    fprintf('  Eb/N0 = %f dB, Simulated BER = %e, Theoretical BER = %e\n', ...
        EbN0_dB, BER(i), theor_BER(i));
end

figure;
semilogy(EbN0_dB_range, BER, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Simulated BER');
hold on;
semilogy(EbN0_dB_range, theor_BER, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical BER (Approx.)');
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-7 1]);
legend('Location', 'southwest');

Ber_rectangular = BER; 
save('Berrectangular.mat', 'Ber_rectangular');

%% 3. Plot the noiseless constellation, both considering only the optimal sampling and the full trajectory
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

%% 4. Plot the constellation at BER target (10^-4), considering only the optimal sampling

BER_target = 1e-4;

% Find the Eb/N0 that results in a BER close to the target from the simulation
% Find the index where simulated BER is closest to the target
[~, target_idx] = min(abs(BER - BER_target));
EbN0_target_dB = EbN0_dB_range(target_idx);

fprintf('Plotting constellation at BER target of %e (approx. Eb/N0 = %f dB)...\n', BER_target, EbN0_target_dB);

EbN0_target_linear = 10^(EbN0_target_dB / 10);
EsN0_target_linear = k * EbN0_target_linear;
sigma_n_squared_target = Ps*SpS./(2*k*EbN0_target_linear);
sigma_n_target = sqrt(sigma_n_squared_target); % Standard deviation

noise_I_plot = sigma_n_target * randn(size(tx_complex_envelope));
noise_Q_plot = sigma_n_target * randn(size(tx_complex_envelope));
rx_complex_envelope_plot =ifft(fft(tx_complex_envelope + (noise_I_plot + 1i * noise_Q_plot)).*H');

% Optimal Sampling
rx_sampled_symbols_plot = rx_complex_envelope_plot(SpS/2:SpS:end);

figure;
plot(real(rx_sampled_symbols_plot), imag(rx_sampled_symbols_plot), 'x', 'MarkerSize', 6, 'Color', [0.85 0.325 0.098]); 
grid on;
xlabel('In-Phase (I)');
ylabel('Quadrature (Q)');
max_val_noisy = max(abs(constellation_points(:))) * 1.5;
xlim([-max_val_noisy max_val_noisy]);
ylim([-max_val_noisy max_val_noisy]);

%% eyediagram
eyediagram(rx_complex_envelope_noiseless(1:10000), 2*SpS, 2*SpS); 