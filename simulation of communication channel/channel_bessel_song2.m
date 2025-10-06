close all
clear all
clc

[data, Fs] = audioread('ABBA - Mamma Mia.mp3'); 
N= length(data);
x_ref_1= data(1:N, 1);
DeltaTime = 1 / Fs;
t1 = (0:DeltaTime:(N-1)*DeltaTime)';

load('channel.mat');
N = length(channel);
x_ch = channel(1:N, 1);
DeltaTime = 1 / Fs;
t = (0:DeltaTime:(N-1)*DeltaTime)';
DeltaFrequency=Fs/N;
f=(-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kHz

%% Receiving signal using Bessel low-pass filter 
fc2 = 12000; % carrier frequency of the second song
x_in = x_ch;
x1_mod = x_ref_1.*cos(2 * pi * fc2 * t1);
x_in(1:length(x1_mod)) = x_in(1:length(x1_mod)) + x1_mod;
x_in = 2*x_in.*cos(2*pi*fc2*t);
x_ref = x_ch.*0;
x_ref(1:length(x_ref_1)) = x_ref_1;

% SIR Optimization (Cutoff Frequency Tuning)
freq_range = 4000:100:8000; % I used frequency range with big steps to reduce processing time but any step can be used.
best_SIR = 0; % intial value for SIR
best_cutoff = 0;
SIR_values = zeros(size(freq_range));
order = 4;                % 4-pole Bessel filter
for k = 1:length(freq_range)
    fc = freq_range(k);
    [b, a] = besself(order, 2*pi*fc);
    [Bz, Az] = bilinear(b, a, Fs);
    x_out = filter(Bz, Az, x_in);
    delay = finddelay(x_ref, x_out);
    x_out = circshift(x_out, -delay);
    e = x_out - x_ref;
    SIR = sum(abs(x_out).^2) / sum(abs(e).^2);
    SIR_values(k) = 10*log10(SIR); % in dB
    if SIR > best_SIR
        best_SIR = SIR;
        best_cutoff = fc;
    end
end

figure(1);
plot(freq_range/1000, SIR_values, '-o');
% title('SIR vs Cutoff Frequency for 4-pole Bessel Filter');
xlabel('Cutoff Frequency (kHz)');
ylabel('SIR (dB)');
grid on;
hold on;
plot(best_cutoff/1000, 10*log10(best_SIR), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(best_cutoff/1000, 10*log10(best_SIR), sprintf('Peak at %.4f kHz', best_cutoff/1000), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
hold off;

fprintf('Optimized Cutoff Frequency: %d Hz\n', best_cutoff);



% tranfer function of optimum filter

[b, a] = besself(order, 2*pi*best_cutoff); % Analog transfer function

[bd, ad] = bilinear(b, a, Fs);   % Convert to Digital using Bilinear Transformation

H_s = tf(b, a);          
H_z = tf(bd, ad, 1/Fs);  

f_range = logspace(0, log10(Fs/2), 1000); % Frequency range in log scale
w = 2 * pi * f_range;

% Calculate Frequency Response (Analog)
[mag, ~] = bode(H_s, w);
mag = squeeze(20*log10(mag)); 

% Convert Frequency Axis to kHz
f_range_kHz = f_range / 1000;
fc_kHz = best_cutoff / 1000;

% Plot Magnitude Response of Analog Transfer Function
figure(2);
semilogx(f_range_kHz, mag, 'b-', 'LineWidth', 2);
hold on;
xline(fc_kHz, '--r', 'Cutoff Frequency', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', ...
    'Color', 'r', 'LineWidth', 1.5);
text(fc_kHz, -7, sprintf('f_c = %.3f kHz', fc_kHz), 'FontSize', 10, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'right');
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');


% comparison plot of filtered and ref signal

fft_original = fftshift(fft(x_ref));

[b, a] = besself(order, 2*pi*best_cutoff);
[Bz, Az] = bilinear(b, a, Fs);
x_out = filter(Bz, Az, x_in);
delay = finddelay(x_ref, x_out);
x_out = circshift(x_out, -delay);
fft_filtered = fftshift(fft(x_out));

figure(3);
plot(f, 10*log10(abs(fft_original).^2), 'b'); hold on;
plot(f, 10*log10(abs(fft_filtered).^2), 'r');
xline(7.27, '--k', 'Bandwidth Limit','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
ylabel('|X(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on;
legend('Original','Filtered (Optimized)');