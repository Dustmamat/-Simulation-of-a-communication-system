clc;
clear;
close all;

nbit = 2;

Ntransmitted = 2^12;


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


fc = 0.3;
H = fftshift(1./(1+1j.*f./fc));

EbN0 = 14.9;
Ps = mean(abs(s.^2))/nbit;

%% Considering noise
    varn = Ps/(2*10^(EbN0/10))*SpS;
    sigma = sqrt(varn);
    n = sigma*randn(Ntransmitted*SpS/nbit,1);
    x = s+n;
    X = fft(x);
    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));

eyediagram(y,2*SpS,2*SpS);

%% Noiseless
    n = 0;
    x = s+n;
    X = fft(x);
    %Filter at the receiver
    Y = X.*H.';
    y = real(ifft((Y)));

eyediagram(y,2*SpS,2*SpS);