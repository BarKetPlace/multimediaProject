function [MFCCcell] = getMFCC(DATA,snr,denoise_flag)
%% Mel filter bank (Tiphanie report section C.2.3 - Filter Bank)

% %Overlap between triangles, percentage of step in mel domain (between 0 and 1)
% Overlap = .5;
 Fs = 16000;
% 
% frameLength_time = 25; %Frame length in ms
% frameLength = frameLength_time/1000*Fs;%Frame length in samples
% DFTlength = frameLength;
% [FilterBank] = MelCepstrumFilterBank(Fs, Overlap, DFTlength);
% 

D=[];

%Load Dictionnary in case of denoising
if denoise_flag
    load Codebooks.mat
    D = Codebooks{1,2}; %We arbitrarly choose a dictionary
%     D=D./(ones(size(D,1),1)*sqrt(sum(D.^2)));
end
   
%load Noise in case of snr != -1
if snr~=-1
    noise_path = '../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';

    [Noise] = audioread([noise_path noise_file]);

end


NbFiles = length(DATA.utt);
%  fprintf('MFCC Extraction:     \n');
for ifile = 1:NbFiles
    fprintf('File %d/%d\n',ifile,NbFiles);
    x = DATA.rawSpeech{1,ifile};
    
%     x= x(DATA.speechframes{ifile});
    if snr~=-1
        %Create noise
        %extract the right noise length
        noise_sig = Noise(1:length(x))';
        %Uniformization of noise
        UVnoise_sig = noise_sig/std(noise_sig);
        UVnoise_sig = UVnoise_sig -mean(UVnoise_sig);
        %future variance of noise depending on SNR
        varn= (var(x)) / (10^(snr/10));
        %Amplification of noise
        n = (varn^(.5))*UVnoise_sig;
        y = x + n;
        y=y(:);
        
        %First we process the clean signal
        [~,Ex,pspectrum] = melfcc(x,[],Fs,[]);
        %Then the noisy signal
        [~,Ey,pspectrum] = melfcc(y,[],Fs,[]);
        
        %Then we use the noise's features to denoise
        [cepstra,aspectrum,pspectrum] = melfcc(y,Ey-Ex,Fs,D);
    else %if snr==-1 (no additive noise)
        y=x;
        [cepstra,aspectrum,pspectrum] = melfcc(y,[],Fs,D);
    end
%     while (m ~= SigLength)
%         yf = y(n:m);
%        
%         MFCC(1:13,iframe) = getFrameMFCC(yf,FilterBank)';
% %         MFCC(1:13,iframe) = cepstra';
%         n = n + frameLength;
%         m = min(SigLength, m+frameLength);
%         iframe=iframe + 1;
%     end
    
    DATA.mfcc{ifile} = cepstra;
end


MFCCcell = DATA.mfcc;
end
