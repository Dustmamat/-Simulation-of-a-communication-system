clear all
close all
clc


%% Analyze 2-signal with disturbance

[data, Fs] = audioread('ABBA - Mamma Mia.mp3'); 
N= length(data);
x_ref_2 = data(1:N, 1);

load('307576_SONG2.mat');
N = length(xtot);
x_2 = xtot(1:N, 1);
duration = N / Fs;

DeltaTime = 1 / Fs;
t = (0:DeltaTime:(N-1)*DeltaTime)';
DeltaFrequency=Fs/N;
f=(-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kHz


%plot of of the signal
figure(1)
plot(t, x_2,'r');
xlabel('time (s)')
ylabel('x(t)')
ylim([-1,1])
grid on


X=fftshift(fft(x_2)); % DTF of the signal

figure(2)
plot(f, 10*log10(abs(X).^2), 'b')
ylabel('|X(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on


%% Bessel low-pass filter for 1-signal
x_in = x_2;
x_ref = x_ref_2;
% SIR Optimization (Cutoff Frequency Tuning)
freq_range = 1000:100:10000; % I used frequency range with step size of 100 to reduce processing time but any step can be used.
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

figure(3);
plot(freq_range/1000, SIR_values, '-o');
% title('SIR vs Cutoff Frequency for 4-pole Bessel Filter');
xlabel('Cutoff Frequency (kHz)');
ylabel('SIR (dB)');
grid on;
hold on;
plot(best_cutoff/1000, 10*log10(best_SIR), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(best_cutoff/1000, 10*log10(best_SIR), sprintf('Peak at %.1f kHz', best_cutoff/1000), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
hold off;

fprintf('Optimized Cutoff Frequency: %d Hz\n', best_cutoff);
fprintf('Best SIR value: %d dB\n', 10*log10(best_SIR));


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
figure(4);
semilogx(f_range_kHz, mag, 'b-', 'LineWidth', 2);
hold on;
xline(fc_kHz, '--r', 'Cutoff Frequency', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', ...
    'Color', 'r', 'LineWidth', 1.5);
text(fc_kHz, -7, sprintf('f_c = %.1f kHz', fc_kHz), 'FontSize', 10, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'right');
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');



% Plot Magnitude Response of Digital Transfer Function
[H, fre] = freqz(bd, ad, 1024, Fs);
f_kHz = fre / 1000;

figure(5);
semilogx(f_kHz, 20*log10(abs(H)), 'b-', 'LineWidth', 2);
hold on;
xline(fc_kHz, '--r', 'Cutoff Frequency', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', ...
    'Color', 'r', 'LineWidth', 1.5);
text(fc_kHz, -10, sprintf('f_c = %.1f kHz', fc_kHz), 'FontSize', 10, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'right');
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
set(gca, 'XScale', 'log'); 


% comparison plot of filtered and ref signal

fft_original = fftshift(fft(x_ref));

[b, a] = besself(order, 2*pi*best_cutoff);
[Bz, Az] = bilinear(b, a, Fs);
x_out = filter(Bz, Az, x_in);
delay = finddelay(x_ref, x_out);
x_out = circshift(x_out, -delay);
fft_filtered = fftshift(fft(x_out));

figure(6);
plot(f, 10*log10(abs(fft_original).^2), 'b'); hold on;
plot(f, 10*log10(abs(fft_filtered).^2), 'r');
ylabel('|X(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on;
legend('Original','Filtered (Optimized)');