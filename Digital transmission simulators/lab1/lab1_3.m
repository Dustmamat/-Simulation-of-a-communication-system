clear all
close all
clc


% Load audio data
[data, Fs] = audioread('Time_10.m4a'); 
N= length(data);
x = data(1:N, 1)'; 
xmax = max(x);
xmin = min(x);
x_mean = mean(x); 
x = x-x_mean; % shift x so that mean value is zero

% uniform quantization 
n = 5; 
M = 2^n; % number of quantization levels
delta = (xmax - xmin) / M;
partition = xmin + delta * (1:M-1); 
codebook = xmin + delta * (0.5 : 1 : M - 0.5); % midpoints of the quantization intervals
[indexes, xq] = quantiz(x, partition, codebook);

% Huffman coding
symbols = 0:M-1;   % the integer numbers [0;M-1] pointing to intervals for each element of Signal
occurances = zeros(1, length(symbols));
for k = 1:length(symbols)
    occurances(k) = sum(indexes == symbols(k));
end

valid = occurances > 0; % Remove zero-probability symbols
symbols = symbols(valid);
occurances = occurances(valid);
p = occurances./sum(occurances);
H = sum(p.*log2(1./p)) % entropy of the source 
Dictionary = huffmandict(symbols, p);
L_fixed = log2(M);
avg_L_huff = 0;
for i = 1:length(Dictionary)
    avg_L_huff = avg_L_huff + p(i) * length(Dictionary{i,2});
end

coding_gain = L_fixed/avg_L_huff;
fprintf('Coding gain: %.2f\n', coding_gain);


figure;
bar(symbols, p);
xlabel('Quantization Level (Symbol)');
ylabel('Probability');
grid on;


