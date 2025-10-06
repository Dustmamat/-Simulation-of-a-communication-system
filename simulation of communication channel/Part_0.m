close all
clear all
clc


%% Analyze 1-reference signal
[data, Fs] = audioread('John Lennon - Imagine.mp3'); 
N= length(data);
x_ref_1 = data(1:N, 1);
DeltaTime = 1 / Fs;
t = (0:DeltaTime:(N-1)*DeltaTime)';
DeltaFrequency=Fs/N;
f = (-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kHz

% Bandwidth defined as frequency range where 95% of energy is contained.
X = fft(x_ref_1);
X = abs(X(1:N/2+1)); % Take positive frequencies only
f_pos = linspace(0, Fs/2, length(X));
S = X.^2; % energy spectral density
E_total = sum(S);
cum_E = cumsum(S); % cumulative energy 
index_95 = find(cum_E >= 0.95 * E_total, 1, 'first');
BW_95 = f_pos(index_95);
fprintf('95%% Energy Bandwidth for 1-ref signal: %.2f Hz\n', BW_95);

% time domain plot of the 1-ref signal
figure(1)
plot(t, x_ref_1,'r');
xlabel('time (s)')
ylabel('x(t)')
ylim([-1,1])
% title('1-ref signal in time domain')
grid on

% signal in frequency domain DTF
X=fftshift(fft(x_ref_1));
figure(2)
plot(f, 10*log10(abs(X)), 'b')
% title('1-ref signal in frequency domain')
ylabel('|X(k)| (dB)')
xlabel('frequency (kHz)')
grid on

%% Analyze 2-ref signal

[data, Fs] = audioread('ABBA - Mamma Mia.mp3'); 
N= length(data);
x_ref_2 = data(1:N, 1);

DeltaTime = 1 / Fs;
t = (0:DeltaTime:(N-1)*DeltaTime)';
DeltaFrequency=Fs/N;
f=(-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kH

% Bandwidth defined as frequency range where 95% of energy is contained.
X = fft(x_ref_2);
X = abs(X(1:N/2+1)); % Take positive frequencies only
f_pos = linspace(0, Fs/2, length(X));
S = X.^2; % energy spectral density
E_total = sum(S);
cum_E = cumsum(S); % cumulative energy 
index_95 = find(cum_E >= 0.95 * E_total, 1, 'first');
BW_95 = f_pos(index_95);
fprintf('95%% Energy Bandwidth for 2-ref signal: %.2f Hz\n', BW_95);

% time domain plot of the 1-ref signal
figure(3)
plot(t, x_ref_2,'r');
xlabel('time (s)')
ylabel('x(t)')
ylim([-1,1])
% title('2-ref signal in time domain')
grid on

% signal in frequency domain
X=fftshift(fft(x_ref_2));
figure(4)
plot(f, 10*log10(abs(X)), 'b')
% title('2-ref signal in frequency domain')
ylabel('|X(k)| (dB)')
xlabel('frequency (kHz)')
grid on