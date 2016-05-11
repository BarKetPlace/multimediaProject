function makenoise(SNR)
%function makenoise(SNR)
%INPUT: SNR you want for the noisy signal
% The function loads dataTest.mat, corrupts it with white noise and save the
% result into dataTestNoisySNRdB.mat
%OUTPUT: none


fprintf(['Creating dataTestNoisy' num2str(SNR) 'dB.mat...']);
load dataTest.mat

noise_path = '/home/antoine/Documents/multimediaProject/TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
[Noise, Fs] = audioread([noise_path noise_file]);
% Noise=Noise-mean(Noise);


% soundsc(Noise,Fs);

for ifile = 1:length(DATA.rawSpeech)
    x = DATA.rawSpeech{1,ifile};
    %Conditioning of x
    x=x-mean(x);
    x=x/max(abs(x));
    %extract the right noise length
    
    noise_sig = Noise(1:length(x))';
    UVnoise_sig = noise_sig/std(noise_sig);
    varn= (var(x)) / (10^(SNR/10));
    noise_sig = (varn^(.5))*UVnoise_sig;
    noise_sig = noise_sig-mean(noise_sig);
%     noise_sig=noise_sig'/sqrt(varn);
    DATA.rawSpeech{1,ifile} = x + noise_sig;
   
end

save(['dataTestNoisy' num2str(SNR) 'dB.mat'],'DATA');
fprintf('done.\n');
end
