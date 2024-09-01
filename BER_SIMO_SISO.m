%% this code is to perform BER analysis of siso with varaiable M
%% 19- august  -2024
tic
clc;
clear all;
close all;

M=4; % modulation order 

N=100000; % Number of data( no of samples)

Nr=4;% number of receiver antennas

input_data = randi([0 M-1], 1, N); % generate random number of integer 0 to M-1 and dimension is 1 row and N columns

bit_input=de2bi(input_data,log2(M),'left-msb');

x=pskmod(input_data,M);

%x=qammod(input_data,M,'UnitAveragePower',true);% consider unit avg power as 1
%x=qammod(input_data,M);

n=sqrt(1/2)*(randn(Nr,N)+1i*randn(Nr,N));% noise term
h=sqrt(1/2)*(randn(Nr,N)+1i*randn(Nr,N));% channel 
SNR=0:2:20;
%SNR=0;
snr=10.^(SNR/10);

for jj=1:1:length(SNR)

y=sqrt(snr(jj))*h.*x+n;% received signal with noise

%y=sqrt(snr(jj))*h.*x; % received signal without noise

x1=pskmod([0:1:M-1],M);% generating a data base or lookup table


for  j1=1:1:N

    for j2=1:1:M
        d(j1,j2) =norm(y(:,j1)-(sqrt(snr(jj))*h(:,j1).*x1(j2))).^2;
    end

end
for j1=1:1:N
    [y2(j1),ind(j1)]=min(d(j1,:));
end
deteceted_symbol=ind-1;
bit_output=de2bi(deteceted_symbol,log2(M),'left-msb');
[berr(jj),bratio(jj)]=biterr(bit_output,bit_input);
[serr(jj),sratio(jj)]=symerr(deteceted_symbol,input_data);
end
semilogy(SNR,sratio,'m-',SNR,bratio,'r-');
xlabel('SNR(db)')
ylabel('BER')
grid on
legend('SER','BER')
toc