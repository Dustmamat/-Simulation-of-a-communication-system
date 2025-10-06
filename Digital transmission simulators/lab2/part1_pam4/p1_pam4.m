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

H = fftshift(sinc(f).*exp(-1j*2*pi*f/2*nbit));

EbN0 = 1:0.5:22; %in dB
BER = zeros(length(EbN0),1);
nError = zeros(length(EbN0),1);

Ps = mean(abs(s.^2))/nbit;
%%
for i=1:length(EbN0)
    varn = Ps/(2*10^(EbN0(i)/10))*SpS/nbit;
    sigma = sqrt(varn);
    n = sigma*randn(Ntransmitted*SpS/nbit,1);
    %Enter the channel
    x = s+n;
    X = fft(x);

    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));
    ysampled=y(SpS:SpS:end);
    for j=1:length(ysampled)
        if(ysampled(j)<=-2)
            ysampled(j) = -3;
        elseif((-2<ysampled(j)) && (ysampled(j)<0))
                ysampled(j) = -1;
        elseif((0<=ysampled(j)) && (ysampled(j)<2))
            ysampled(j) = 1;
        else
            ysampled(j) = 3;
        end
    end
    nError(i) = length(find(ysampled-s(SpS:SpS:end)));
    BER(i) = nError(i)/Ntransmitted;
    display("Running"+num2str(i/length(EbN0)*100)+"% completed");

end

%%
M = 4;
figure, hold on;
% Theory
BER_theory =(M-1)/(M*log2(M))*erfc(sqrt(3*log2(M)/(M^2-1)*10.^(EbN0/10)));
%BER_theory = 3/8*erfc(sqrt(2/5*10.^(EbN0/10)))';

plot(EbN0,BER_theory,Linewidth=1.5);

plot(EbN0,BER,'x',LineWidth=2);
yscale log;
grid on;
title("PAM-4 BER against EbN0");
xlabel("$\frac{E_b}{N_0}$ (dB)",interpreter="latex");
ylabel("BER",interpreter='latex');
xlim([0 15.5])
legend("Theory","Simulation");