
close all;
clear all;
clc;

% Load BER data for all 8-QAM types
load('Berrectangular.mat', 'Ber_rectangular');
load('Bercircular.mat', 'Ber_circular');
load('Bersp.mat', 'Ber_sp');
load('BerX.mat', 'Ber_X');
load('Bersquare.mat', 'Ber_square');
load('Berstar.mat', 'Ber_star');

EbN0_dB_range = 2:1:15; % Eb/N0 range in dB

figure;
% Generate 6 distinct colors for the curves
colors = lines(6); 

% Plot each BER curve with a unique color 
semilogy(EbN0_dB_range, Ber_X, 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', '8X - 8QAM'); hold on;
semilogy(EbN0_dB_range, Ber_star, 'Color', colors(2,:), 'LineWidth', 1.5, 'DisplayName', 'Star 8-QAM');
semilogy(EbN0_dB_range, Ber_square, 'Color', colors(3,:), 'LineWidth', 1.5, 'DisplayName', 'Square 8-QAM');
semilogy(EbN0_dB_range, Ber_rectangular, 'Color', colors(4,:), 'LineWidth', 1.5, 'DisplayName', 'Rectangular 8-QAM');
semilogy(EbN0_dB_range, Ber_sp, 'Color', colors(5,:), 'LineWidth', 1.5, 'DisplayName', 'SP 8-QAM');
semilogy(EbN0_dB_range, Ber_circular, 'Color', colors(6,:), 'LineWidth', 1.5, 'DisplayName', 'Circular 8-QAM (8-PSK)');

grid on; 
xlabel('Eb/N0 (dB)'); 
ylabel('Bit Error Rate (BER)'); 
legend('Location', 'southwest');
ylim([1e-6 1]); 
xlim([min(EbN0_dB_range) max(EbN0_dB_range)]);
hold off; 