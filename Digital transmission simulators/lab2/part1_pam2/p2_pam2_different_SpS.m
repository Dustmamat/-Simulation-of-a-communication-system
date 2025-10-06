clc;
clear;
close all;

nbit = 1;
Ntransmitted = 2^20;
B = randi([0 1],[Ntransmitted,1]);

SpS = 8;

f =(-SpS*Ntransmitted/2:SpS*Ntransmitted/2-1)/Ntransmitted;
%Define the Antipodal PAM-2 Alphabet
alfa = [-1, 1];

%Now to formulate the electrical signal
s = zeros(SpS*Ntransmitted,1); %if not defined beforehand can slow down the program a lot!!


for i=1:Ntransmitted
    s((i-1)*SpS+1:i*SpS) = ones(1,SpS)*alfa(B(i)+1);
end

fc = 0.4; %choose the optimal fc
H = fftshift(1./(1+1j.*f./fc));

EbN0 = 8:0.1:20;
BER = zeros(length(EbN0),length(fc));
nError = zeros(length(EbN0),length(fc));
%%
Ps = 1;
sampling_instant = 1:SpS;
for i=1:length(EbN0)
    varn = Ps/(2*10^(EbN0(i)/10))*SpS;
    sigma = sqrt(varn);
    n = sigma*randn(Ntransmitted*SpS,1);
    x = s+n;
    X = fft(x);
    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));
    for jj=1:length(sampling_instant)
    ysampled=y(sampling_instant(jj):SpS:end);
    for k=1:Ntransmitted
        if(ysampled(k)>0)
            ysampled(k) = 1;
        else
            ysampled(k)=-1;
        end
    end
    nError(i,jj) = length(find(ysampled-s(sampling_instant(jj):SpS:end)));
    BER(i,jj) = nError(i,jj)/Ntransmitted;

    end
    display("Running"+num2str(i/length(EbN0)*100)+"% completed");
end
%%
figure, hold on;
% Theory
BER_theory = 0.5 .* erfc(sqrt(10.^(EbN0/10)))';
plot(EbN0, BER_theory, 'LineWidth', 1.5);

% Initialize legend entries
legend_entries = ["Theory-Matched Filter"];

% Plot simulated BER for each fc
for i = 1:length(sampling_instant)
    plot(EbN0, BER(:,i), 'LineWidth', 2);
    legend_entries(end+1) = "Instant = " + num2str(sampling_instant(i));
end

% Configure plot
set(gca, 'YScale', 'log');
grid on;
title("PAM-2 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)", 'Interpreter', "latex");
ylabel("BER", 'Interpreter', 'latex');
ylim([1.52272e-10,1]);
xlim([8,11.55446583990919]);

legend(legend_entries, 'Location', 'best');



%% Now to plot fc against EbN0 required for BER TARGET = 1e-4
EbN0_req = zeros(length(fc),1);

for i=1:length(fc)
    v = find(BER(:,i)<=1e-4);
    if(isempty(v))
        EbN0_req(i) = Inf;
    else
        EbN0_req(i) = EbN0(v(1));
    end
end
figure;
plot(fc,EbN0_req);