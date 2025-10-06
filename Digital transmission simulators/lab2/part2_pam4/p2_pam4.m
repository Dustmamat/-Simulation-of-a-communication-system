clc;
clear;
close all;

nbit = 2;

Ntransmitted = 2^20;


B = randi([0 1],[Ntransmitted,1]);

SpS = 8;
alfa = [-3 -1 3 1]; %the reason why the 3 and 1 are swapped is to enable gray coding
pam4_dict = dictionary([0 1 2 3],alfa);
symbol_index = @(bit1,bit2) 2*bit1+bit2;

%if not defined beforehand can slow down the program a lot!!

%Formulate the electrical signal
s = zeros(SpS*Ntransmitted/nbit,1);
for i=1:nbit:Ntransmitted
    s((i-1)/nbit*SpS+1:((i+1*(nbit-1))/nbit*SpS)) = ones(1,SpS)*pam4_dict(symbol_index(B(i),B(i+1)));

end
f =(-SpS*Ntransmitted/(2*nbit):SpS*Ntransmitted/(2*nbit)-1)/Ntransmitted;


fc = 0.1:0.1:1.5;
fc = fc';
H = fftshift(1./(1+1j.*f./fc));

EbN0 = 12:0.1:20;
BER = zeros(length(EbN0),length(fc));
nError = zeros(length(EbN0),length(fc));
Ps = mean(abs(s.^2))/nbit;

%%
for i=1:length(EbN0)
    varn = Ps/(2*10^(EbN0(i)/10))*SpS;
    sigma = sqrt(varn);
    n = sigma*randn(Ntransmitted*SpS/nbit,1);
    x = s+n;
    X = fft(x);
    for jj =1:length(fc)
        H = fftshift(1./(1+1j*f/fc(jj)));
    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));
    ysampled=y(SpS:SpS:end);
    for k=1:length(ysampled)
        if(ysampled(k)<=-2)
            ysampled(k) = -3;
        elseif((-2<ysampled(k)) && (ysampled(k)<0))
                ysampled(k) = -1;
        elseif((0<=ysampled(k)) && (ysampled(k)<2))
            ysampled(k) = 1;
        else
            ysampled(k) = 3;
        end
    end
    nError(i,jj) = length(find(ysampled-s(SpS:SpS:end)));
    BER(i,jj) = nError(i,jj)/Ntransmitted;

    end
    display("Running"+num2str(i/length(EbN0)*100)+"% completed");
end
%%
figure, hold on;
% Theory
M = 4;
BER_theory =(M-1)/(M*log2(M))*erfc(sqrt(3*log2(M)/(M^2-1)*10.^(EbN0/10)));
plot(EbN0, BER_theory, 'LineWidth', 1.5);

% Initialize legend entries
legend_entries = ["Matched Filter"];

% Plot simulated BER for each fc
for i = 1:length(fc)/2
    plot(EbN0, BER(:,i), 'LineWidth', 2);
    legend_entries(end+1) = "fc = " + num2str(fc(i));
end

% Configure plot
set(gca, 'YScale', 'log');
grid on;
title("PAM-2 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)", 'Interpreter', "latex");
ylabel("BER", 'Interpreter', 'latex');
ylim([1e-7,1e0]);
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
plot(fc,EbN0_req,Linewidth=2);
title("PAM-4 EbN0 at BER=10^{-4} against Filter -3dB Bandwidth");
xlabel("f_c");
ylabel("EbN0");
grid on;