close all
clear all
clc

% Eyediagram of optimal sigle pole filter in noiseless and at EbN0 when BER = 1e-4;

%% System parameters
SpS = 8;                  % Samples per symbol
Tb = 1;                   % Bit duration
Rs = 1/Tb;                % Symbol rate
Nbits = 2^14;             % Number of bits to simulate
n_pam2 = 1;               % Bits per symbol for PAM-2
N_pam2 = Nbits / n_pam2;  % Calculate number of PAM-2 symbols
EbN0_dB = 8:0.1:13;       % for less computation time we chose this range 
EbN0 = 10.^(EbN0_dB./10); % Linear EbN0
BER_target = 1e-4; 
best_fc = 0.4;
Fs = SpS/Tb;              % sampling frequency
deltaFrequency = Fs/(SpS*Nbits);
f = (-Fs/2:deltaFrequency:Fs/2-deltaFrequency);

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
BER_values = zeros(length(EbN0), 1);

for k = 1:length(EbN0)
    noise = sigma(k) * randn(length(tx), 1);
  
    rx_n = tx + noise; % Add noise to transmitted signal
   
    Rx_n = fft(rx_n);    
    H = 1./(1+1i.*f./best_fc); 
    H = fftshift(H);
    Y_n = Rx_n .* H.';
    y_n = ifft(Y_n);
    
    % Sample at optimal points
    rx_samples = y_n(sample_points);
    
    % Make binary decisions (threshold at 0)
    rx_bits = real(rx_samples) > 0;
    
    BER_values(k) = mean(rx_bits ~= bits_pam2);
    
    display("Running"+num2str(k/length(EbN0)*100)+"% completed");
end

%% eyediagram at EbN0 when BER = 1e-4;
[~, idx] = min(abs(BER_values(:) - BER_target));
noise = sigma(idx) * randn(length(tx), 1);
rx_n = tx + noise;
Rx_n = fft(rx_n);    
H = 1./(1+1i.*f./best_fc); 
H = fftshift(H);
Y_n = Rx_n .* H.';
y_n = real(ifft(Y_n));

eyediagram(y_n, 2*SpS, 2*SpS)

% eyediagram in noiseless case 

rx = tx;
Rx = fft(rx);    
H = 1./(1+1i.*f./best_fc); 
H = fftshift(H);
Y = Rx .* H.';
y = real(ifft(Y));

eyediagram(y, 2*SpS, 2*SpS)