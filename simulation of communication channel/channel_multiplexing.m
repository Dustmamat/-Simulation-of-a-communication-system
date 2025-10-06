close all
clear all
clc

%% Analyze channel 
[data, Fs] = audioread('John Lennon - Imagine.mp3'); 
N= length(data);
x_ref_1= data(1:N, 1);
DeltaTime = 1 / Fs;
t1 = (0:DeltaTime:(N-1)*DeltaTime)';

[data, Fs] = audioread('ABBA - Mamma Mia.mp3'); 
N= length(data);
x_ref_2 = data(1:N, 1);
DeltaTime = 1 / Fs;
t2 = (0:DeltaTime:(N-1)*DeltaTime)';

load('channel.mat');
N = length(channel);
x_ch = channel(1:N, 1);

DeltaTime = 1 / Fs;
t = (0:DeltaTime:(N-1)*DeltaTime)';
DeltaFrequency=Fs/N;
f=(-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kHz

% Bandwidth defined as frequency range where 95% of energy is contained.
X = fft(x_ch);
X = abs(X(1:N/2+1)); % Take positive frequencies only
f_pos = linspace(0, Fs/2, length(X));
S = X.^2; % energy spectral density
E_total = sum(S);
cum_E = cumsum(S); % cumulative energy 
index_95 = find(cum_E >= 0.95 * E_total, 1, 'first');
BW_95 = f_pos(index_95);
fprintf('95%% Energy Bandwidth for channel signal: %.2f Hz\n', BW_95);


% %plot of the signal
% figure(1)
% plot(t, x_ch,'r');
% xlabel('time (s)')
% ylabel('x(t)')
% ylim([-3,3])
% grid on


X=fftshift(fft(x_ch)); % DTF of the signal

figure(2)
plot(f, 10*log10(abs(X).^2), 'b')
ylabel('|X(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on

% transmission over the channel without modulation

x_in = x_ch; 
x_in(1:length(x_ref_1)) = x_in(1:length(x_ref_1)) + x_ref_1;
x_out = x_in; 
e = x_ch;
SIR = sum(abs(x_out).^2) / sum(abs(e).^2);
fprintf('SIR of the transmitted song: %.2f dB\n', SIR);

%% Multiplexing two songs

fc1 = 2000;
fc2 = 12000;
x1_mod = x_ref_1.*cos(2 * pi * fc1 * t1);
x2_mod = x_ref_2.*cos(2 * pi * fc2 * t2);
x_in = x_ch;
x_in(1:length(x1_mod)) = x_in(1:length(x1_mod)) + x1_mod;
x_in(1:length(x2_mod)) = x_in(1:length(x2_mod)) + x2_mod; % x_in = x_ch + x1_mod+ x2_mod
transmitted_signal = x_in;

Y=fftshift(fft(transmitted_signal)); % DTF of the signal

figure(3)
plot(f, 10*log10(abs(Y).^2), 'b')
xline(160/1000, '--k', 'Bandwidth of baseband signal','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
xline((fc1-970)/1000, '--k', 'Bandwidth of first song','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
xline((fc1+970)/1000, '--k', 'Bandwidth of first song','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
xline((fc2-7270)/1000, '--k', 'Bandwidth of secong sono','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
ylabel('|Y(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on
