clear all
close all
clc

[data, Fs] = audioread('John Lennon - Imagine.mp3'); 
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

fc1 = 10000; % carrier frequency of the first song
x_in = x_ch;
x1_mod = x_ref_1.*cos(2 * pi * fc1 * t1);
x_in(1:length(x1_mod)) = x_in(1:length(x1_mod)) + x1_mod;
x_in = 2*x_in.*cos(2*pi*fc1*t);
x_ref = x_ch.*0;
x_ref(1:length(x_ref_1)) = x_ref_1;

fre_range = 1000:50:10000;
best_SIR = 0; 
best_f = 0;
k = 0:N-1;
SIR_values = zeros(size(fre_range));
R_p = 0.8;
for k = 1:length(fre_range)
    phi = 2*pi*fre_range(k)/Fs;
    z_1 = 1*exp(1j*2*pi*fc1/Fs);
    p_1 = R_p*exp(1j*phi);
    a = [1, -(p_1 + conj(p_1)), p_1*conj(p_1)];
    b = [1, -(z_1+conj(z_1)), z_1*conj(z_1)];
    [H, W] = freqz(b, a, 1024, Fs);
    H_ref = interp1(W, abs(H), 0, 'linear'); % to renormalize the signal we used zero frequency gain
    x_out = filter(b, a, x_in)/H_ref;

    % here I tried to calculte using DFT but my laptop doesnt allow me to
    % do it because array size is limited. Instead it shows the foollow
    % Requested array exceeds the maximum possible variable size.
    % Error in untitled4 (line 22)
    %     X_out = fft(x_in).*H_k;
    % 
    % Related documentation

    % z = exp(-1j * 2 * pi .* k / N); % obtaining DFT from z-transform
    % H_k = (z-z_1).*(z-conj(z_1))./((z-p_1).*(z-conj(p_1)));

    % X_out = fft(x_in).*H_k;
    % x_out = ifft(X_out);

    e = x_ref - x_out; % error
    SIR = sum(abs(x_out).^2) / sum(abs(e).^2);
    SIR_values(k) = 10*log10(SIR); % in dB
    if SIR > best_SIR
        best_SIR = SIR;
        best_f = fre_range(k);
    end
end


figure(1);
plot(fre_range/1000, SIR_values, '-o');
xlabel('Cutoff Frequency (kHz)');
ylabel('SIR (dB)');
grid on;
hold on;
plot(best_f/1000, 10*log10(best_SIR), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(best_f/1000, 10*log10(best_SIR), sprintf('Peak at %.3f kHz', best_f/1000), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7, 'FontWeight', 'bold');
hold off;

fprintf('Optimized zero frequency: %d Hz\n', best_f);
fprintf('Best SIR value: %d dB\n', 10*log10(best_SIR));

% Transfer function of the optimum filter in Z domain

z1 = 1*exp(1j*2*pi*fc1/Fs);
p1 = R_p*exp(1j*2*pi*best_f/Fs);
z2 = conj(z1);
p2 = conj(p1);
z_step=0.05;
z_limit=1.5;
Re_z_v=[-z_limit:z_step:z_limit];
Im_z_v=[-z_limit:z_step:z_limit];

[Re_z,Im_z]=meshgrid(Re_z_v,Im_z_v);
Z=Re_z+1j*Im_z;

H_z=(Z-z1).*(Z-z2)./(Z-p1)./(Z-p2);

freq=[-0.5:0.001:0.5]*Fs;
z=exp(1j*2*pi*freq/Fs);
H_f=(z-z1).*(z-z2)./(z-p1)./(z-p2);

figure(2);
mesh(Re_z,Im_z,abs(H_z));
grid on;
xlabel('Re(z)');
ylabel('Im(z)');
zlabel('| H(z) |');
hold on;
plot3(real(z),imag(z),abs(H_f),'r-');

figure(3);
xlabel('frequency [kHz]');
ylabel('| H(f) |  [dB]');
hold on;
plot(freq/1000,20*log10(abs(H_f)),'r-');
grid on;

% comparison plot of filtered and ref signal

fft_original = fftshift(fft(x_ref));
z_1 = 1*exp(1j*2*pi*fc1/Fs);
p_1 = R_p*exp(1j*2*pi*best_f/Fs);
a = [1, -(p1 + conj(p1)), p1*conj(p1)];
b = [1, -(z1+conj(z1)), z1*conj(z1)];
[H, W] = freqz(b, a, 1024, Fs);
H_ref = interp1(W, abs(H), 0, 'linear');
x_out = filter(b, a, x_in);
fft_filtered = fftshift(fft(x_out))/H_ref;

DeltaFrequency=Fs/N;
f=(-Fs/2:DeltaFrequency:Fs/2-DeltaFrequency)/1000; % in kHz

figure(4);
plot(f, 10*log10(abs(fft_original).^2), 'b'); hold on;
plot(f, 10*log10(abs(fft_filtered).^2), 'r');
xline(0.97, '--k', 'Bandwidth Limit','LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'FontSize', 5);
ylabel('|X(k)|^2 (dB)')
xlabel('frequency (kHz)')
grid on;
legend('Original','Filtered (Optimized)');


% sound(x_out, Fs)