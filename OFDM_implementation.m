clc; clearvars; close all;

ntap=3;
pdf = [0.5 0.2 0.3];%% power delay profile
pdf = pdf./sqrt(pdf);
N=256;
M = 4;
Cp = N/4;
mcp = zeros(1,1000);
snr = 0:5:25;
k = log2(M);
bitstream = randi([0 1],512,1);
Nbits = length(bitstream);
Pbcp = zeros(1,length(snr));
for i = 1:length(snr)
     for j = 1:1000
         % 3 tap frequency selective channel
        h=sqrt(0.5)*sqrt(1/ntap)*(randn(1,ntap)+ 1j*randn(1,ntap)).*pdf;
        % FFT of channel coefficients
        g = fft(h,N);
        data = reshape(bitstream,Nbits/k,k);
        data_dec = bi2de(data);
        qm = qammod((data_dec),M,'gray');
        ifft_data = ifft(qm,N);
        data_cp = [ifft_data(N-Cp+1:end); ifft_data];%% CP Added
        %%channel1
        rxcp = conv(data_cp,h);
        outcp = awgn(rxcp,snr(i),'measured');
        out_rcp = outcp(Cp+1:end-ntap+1); %% CP Removed
        fft_rcp = fft(out_rcp,N)./g.'; %%Equalization
        q_demodcp = qamdemod(fft_rcp,M); %%QAM Demodulation
        r_bitstreamcp = de2bi(q_demodcp);
        r_bitstreamcp = reshape(r_bitstreamcp,Nbits,1);
        mcp(j) = sum(xor(bitstream,r_bitstreamcp))/Nbits;
     end
    Pbcp(i) = mean(mcp);
end

 figure(1);
 semilogy(snr,Pbcp,'--o');
 xlabel('SNR (dB)');
 ylabel('BER');
 grid on;
 title('OFDM vs BER Comparison for 4-QAM Modulated Signal');
 figure(2);
 plot(abs(g));
 grid on;
 title('Channel Frequency Response for N=256 Points');
 xlim([0 256]);



