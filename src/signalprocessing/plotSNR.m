function [snr_denoise, snr_mel_energy ] = plotSNR(ifig,Ex,Exhat,En)
% function [snr_denoise, snr_mel_energy ] = plotSNR(ifig,Ex,Exhat,En)

snr_mel_energy= 10*log10(sum(Ex.^2)./sum(En.^2));
% snr_mel_energy_model=10*log10(sum(Ex.^2)./sum(En_model.^2));
snr_denoise = 10*log10(sum(Ex.^2)./sum((Ex-Exhat).^2));


figure(ifig); clf; ifig=ifig+1;
plot(snr_denoise); hold on;
plot(snr_mel_energy); hold on;
legend('||Ex||_2/||Ex-Exhat||_2 in dB','||Ex||_2/||En||_2 in dB');
xlabel('Frame number'); ylabel('SNR in dB');
title('SNR compare in Mel energy space');

end

