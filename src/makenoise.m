function makenoise(SNR)

load dataTest.mat

noise_path = '/home/antoine/Documents/multimediaProject/TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
[Noise, Fs] = audioread([noise_path noise_file]);
% soundsc(Noise,Fs);

for ifile = 1:length(DATA.rawSpeech)
    x = DATA.rawSpeech{1,ifile};
    noise_sig = Noise(1:length(x));
%     varn = var(noise_sig);
%     noise_sig=noise_sig'/sqrt(varn);
    DATA.rawSpeech{1,ifile} = 10.^(SNR/10).*x/std(x)+noise_sig'/std(noise_sig);
end

save(['dataTestNoisy' num2str(SNR) 'dB.mat'],'DATA');
end
