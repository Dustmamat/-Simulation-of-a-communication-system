close all
clear all
clc

%% System parameters
SpS = 8;                  % Samples per symbol
Tb = 1;                   % Bit duration
Rs = 1/Tb;                % Symbol rate
Nbits = 2^20;             % Number of bits to simulate
n_pam2 = 1;               % Bits per symbol for PAM-2
N_pam2 = Nbits / n_pam2;  % Calculate number of PAM-2 symbols
EbN0_dB = 8:0.1:13;       % for less computation time we chose this range 
EbN0 = 10.^(EbN0_dB./10); % Linear EbN0
BER_target = 1e-4; 

Fs = SpS/Tb;              % sampling frequency
deltaFrequency = Fs/(SpS*Nbits);
f = (-Fs/2:deltaFrequency:Fs/2-deltaFrequency);


% for less computation time we chose small range ( any range can be)
% relying on theory 
fc = (0.1:0.05:1)'; 

%% PAM-2 Signal Generation
% Generate random bits
bits_pam2 = randi([0 1], Nbits, 1);

% Symbol mapping for PAM-2: 0 → -1, 1 → +1 (antipodal)
symbols_pam2 = 2 * bits_pam2 - 1;

% Generate NRZ pulse
pulse = ones(1, SpS);

% Create transmitted signal 
tx = zeros(N_pam2 * SpS, 1);
for i = 1:N_pam2
    tx((1 + (i-1)*SpS):(i*SpS)) = symbols_pam2(i) * pulse;
end

 
Ps = mean(abs(tx).^2);          % signal power
sigma = sqrt(Ps*SpS./(2*EbN0)); % noise std

%% BER Simulation

% sample at optimal sampling instant
start = SpS;
sample_points = start:SpS:length(tx);

EbN0_fre_values = zeros(length(fc), 1);
BER_values = zeros(length(EbN0), length(fc));

for k = 1:length(EbN0)
    noise = sigma(k) * randn(length(tx), 1);
  
    rx = tx + noise; % Add noise to transmitted signal
   
    Rx = fft(rx);    
    for i = 1:length(fc)
        H = 1./(1+1i.*f./fc(i)); 
        H = fftshift(H);
        Y = Rx .* H.';
        y = ifft(Y);
        
        % Sample at optimal points
        rx_samples = y(sample_points);
        
        % Make binary decisions (threshold at 0)
        rx_bits = real(rx_samples) > 0;
        
        BER_values(k, i) = mean(rx_bits ~= bits_pam2);
    end
    display("Running"+num2str(k/length(EbN0)*100)+"% completed");
end

for i = 1:length(fc)
    [~, idx] = min(abs(BER_values(:, i) - BER_target));
    EbN0_fre_values(i) = EbN0_dB(idx);
end
% Plot results
%% Plot results

figure(1), hold on;
% Matched filter
BER_theory = 0.5 .* erfc(sqrt(EbN0))';
plot(EbN0_dB, BER_theory, 'LineWidth', 1.5);

[~, idx_best_fc] = min(EbN0_fre_values);
best_fc = fc(idx_best_fc);
best_EbN0 = EbN0_fre_values(idx_best_fc);
x1 = 10*log10((erfcinv(2*1e-4))^2);
x2 = EbN0_fre_values(idx_best_fc);
penalty = x2-x1;

% Add text label showing penalty
text(mean([x1, x2]), 1.3e-4, ...
     sprintf('Penalty = %.2f dB', penalty), ...
     'HorizontalAlignment', 'center', ...
     'Interpreter', 'latex');

% Initialize legend entries
legend_entries = ["Matched filter"];

% Plot simulated BER for each fc
for i = 1:3:length(fc)-1
    plot(EbN0_dB, BER_values(:,i), 'LineWidth', 2);
    legend_entries(end+1) = "fc = " + num2str(fc(i));
end

% Configure plot
set(gca, 'YScale', 'log');
grid on;
title("PAM-2 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)", 'Interpreter', "latex");
ylabel("BER", 'Interpreter', 'latex');
legend(legend_entries, 'Location', 'best');

% Add horizontal line at BER = 1e-4
yline(1e-4, '--k', ...
    'Interpreter', 'latex', ...
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

figure(2)
plot(fc, EbN0_fre_values, 'LineWidth', 2);
hold on;

% Add marker at the optimal point
plot(best_fc, best_EbN0, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
text(best_fc, best_EbN0 + 0.3, ...
    sprintf('Best f_c = %.2f', best_fc), ...
    'HorizontalAlignment', 'left');

grid on;
xlabel('Filter Cutoff Frequency');
ylabel('Eb/N0 (dB)');

% Print results
fprintf('Best cutoff frequency: %.3f \n', best_fc);
fprintf('Penalty of single pole filter compared to the matched filter: %.2f dB\n', penalty);
fprintf('Minimum Eb/N0 at BER = 1e-4: %.2f dB\n', min(EbN0_fre_values));