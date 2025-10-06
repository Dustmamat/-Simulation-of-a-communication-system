clc;
clear;
close all;

nbit = 1;
Ntransmitted = 2^19;
B = randi([0 1],[Ntransmitted,1]);

SpS = 8;

f =(-SpS*Ntransmitted/2:SpS*Ntransmitted/2-1)/Ntransmitted;
%Define the Antipodal PAM-2 Alphabet
alfa = [-1, 1];

%Now to formulate the electrical signal
s = zeros(SpS*Ntransmitted,1); %if not defined beforehand can slow down the program a lot!!

H = fftshift(sinc(f).*exp(-1j*2*pi*f/2));
for i=1:Ntransmitted
    s((i-1)*SpS+1:i*SpS) = ones(1,SpS)*alfa(B(i)+1);
end


%%

Ps = 1; %%for pam2

EbN0 = 1:0.3:11; %in dB
BER = zeros(length(EbN0),1);
nError = zeros(length(EbN0),1);
for i=1:length(EbN0)
    varn = Ps/(2*10^(EbN0(i)/10))*SpS;
    sigma = sqrt(varn);
    n = sigma*randn(Ntransmitted*SpS,1);
    %Enter the channel
    x = s+n;
    X = fft(x);

    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));
    ysampled=y(SpS:SpS:end);
    for j=1:Ntransmitted
        if(ysampled(j)>0)
            ysampled(j) = 1;
        else
            ysampled(j)=-1;
        end
    end
    nError(i) = length(find(ysampled-s(SpS:SpS:end)));
    BER(i) = nError(i)/Ntransmitted;
    
    display("Running"+num2str(i/length(EbN0)*100)+"% completed");

end
%% Plots
figure, hold on;
BER_theory = 0.5.*erfc(sqrt(10.^(EbN0/10)))';
plot(EbN0,BER_theory,Linewidth=1.5);
plot(EbN0,BER,'x',LineWidth=2);
yscale log;
grid on;
title("PAM-2 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)",interpreter="latex");
ylabel("BER",interpreter='latex');
legend("Theory","Simulation");
