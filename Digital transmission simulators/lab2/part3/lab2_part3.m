clear all
close all
clc


%% Simulation parameters
SpS = 8;                  % Samples per symbol
Tb = 1;
Nbits = 2^20;             % Number of bits to simulate
n_pam2 = 1;               % Bits per symbol for PAM-2
N_pam2 = Nbits / n_pam2;  % Calculate number of PAM-2 symbols
EbN0_dB = 9:0.1:20; 
EbN0 = 10.^(EbN0_dB./10);
penalty = 0.5; % in dB
BER_target = 1e-4;
best_f = 0.4;             
Vth = -1:0.04:1;       

Fs = SpS/Tb;              % sampling frequency
deltaFrequency = Fs/(SpS*Nbits);
f = (-Fs/2:deltaFrequency:Fs/2-deltaFrequency);

%% Generation of signal 
% Generate random bits
bits_pam2 = randi([0 1], Nbits, 1);

% Symbol mapping for PAM-2: 0 → -1, 1 → +1 (antipodal)
symbols_pam2 = 2 * bits_pam2 - 1;

% Generate NRZ pulse
pulse = ones(1, SpS);

% Create transmitted signal (upsampled and pulse-shaped)
tx = zeros(N_pam2 * SpS, 1);
for i = 1:N_pam2
    tx((1 + (i-1)*SpS):(i*SpS)) = symbols_pam2(i) * pulse;
end


Ps = mean(abs(tx).^2);
sigma = sqrt(Ps*SpS./(2*EbN0)); % noise std

%% Compute filter
H = 1./(1+1i.*f./best_f);
H = fftshift(H);
% Sample points
start = SpS;
sample_points = start:SpS:length(tx);

%% BER Simulation 

EbN0_Vth_values = zeros(length(Vth), 1);
BER_values = zeros(length(EbN0), length(Vth));

for k = 1:length(EbN0)
    
    noise = sigma(k)*randn(length(tx), 1); % Generate AWGN for this EbN0 level
    rx = tx + noise;                       
    Rx = fft(rx);
    Y = Rx.*H.';
    y = ifft(Y);  
    
    % Sample at optimal sampling points
    rx_samples = y(sample_points);
    
    % Check all threshold values for this noise level
    for i = 1:length(Vth)
        % Make binary decisions (threshold at Vth(i))
        rx_bits = real(rx_samples) > Vth(i);
        % Compare with original bits and calculate BER
        BER_values(k, i) = mean(rx_bits ~= bits_pam2);
    end
    display("Running"+num2str(k/length(EbN0)*100)+"% completed");
end

%% find the Eb/N0 when Vth = 0
EbN0_Vth_0 = 0;
idx_Vth_0 = 0;
for i = 1:length(Vth)
    if Vth(i)== 0
    [~, idx] = min(abs(BER_values(:,i) - BER_target));
    EbN0_Vth_0 = EbN0_dB(idx);
    idx_Vth_0 = i;
    break
    end
end
%% Find the Eb/N0 that gives closest to target BER for each Vth
for i = 1:length(Vth)
    idx_list = find(BER_values(:, i) <= BER_target);
    if ~isempty(idx_list)
        idx = idx_list(1);  % First Eb/N0 that achieves BER <= target
        EbN0_Vth_values(i) = EbN0_dB(idx);
    else
        EbN0_Vth_values(i) = NaN;  % Or Inf, or some flag to indicate no match
    end
end
% find the Vth imposing 0.5 dB penalty
diffs = abs((EbN0_Vth_values(idx_Vth_0+1:end) - EbN0_Vth_0) - penalty);
min_val = min(diffs);
all_matches = find(diffs == min_val);  % All matching indices
first_match = all_matches(1);       

Vth_max = Vth(idx_Vth_0 + first_match);

%% Plot results

figure(1), hold on;
% Matched filter
m = 31;
BER_theory = 0.5 .* erfc(sqrt(EbN0(1:m)))'; % take few values of EbN0
plot(EbN0_dB(1:m), BER_theory, 'LineWidth', 1.5);


% Initialize legend entries
legend_entries = ["Matched filter"];

% Plot simulated BER for five Vths
indices = 20:2:32;

for idx = indices
    plot(EbN0_dB(1:m), BER_values((1:m), idx), 'LineWidth', 2);
    legend_entries(end+1) = "Vth = " + num2str(Vth(idx));
end

% Configure plot
set(gca, 'YScale', 'log');
grid on;
title("PAM-2 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)", 'Interpreter', "latex");
ylabel("BER", 'Interpreter', 'latex');
legend(legend_entries, 'Location', 'best');
% % Add horizontal line at BER = 1e-4
yline(1e-4, '--k', ...
    'Interpreter', 'latex', ...
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

x1 = EbN0_Vth_values(idx_Vth_0)+penalty;
x2 = EbN0_Vth_values(idx_Vth_0);
xline(x1, '--k', 'HandleVisibility', 'off', 'Interpreter', 'latex');
xline(x2, '--k', 'HandleVisibility', 'off', 'Interpreter', 'latex');


% Add text label showing penalty
text(mean([x1, x2]), 1.5e-5, ...
     sprintf('Penalty = %.2f dB', penalty), ...
     'HorizontalAlignment', 'center', ...
     'Interpreter', 'latex');

figure(2)
plot(Vth, EbN0_Vth_values,'LineWidth', 3);
grid on;
xlabel('Decision Threshold ');
ylabel('Eb/N0 (dB)');
% Add vertical lines at Vth_max
xline(Vth_max, '--r', ...
    ['V_{th\_max} = ', num2str(Vth_max, '%.2f')], ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'Color', 'r', ...
    'HandleVisibility', 'off'); 

xline(-Vth_max, '--r', ...
    ['-V_{th\_max} = ', num2str(-Vth_max, '%.2f')], ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'Color', 'r', ...
    'HandleVisibility', 'off'); 
% Print results
fprintf('Maximum threshold error imposing 0.5 dB penalty: %.3f \n', Vth_max);
