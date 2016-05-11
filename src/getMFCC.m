function [MFCCcell] = getMFCC(DATA,denoise_flag)
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
    D = Codebooks{1,5}; %We arbitrarly choose a dictionary
end

NbFiles = length(DATA.utt);
%  fprintf('MFCC Extraction:     \n');
for ifile = 1:NbFiles
%     fprintf('file %d/NbFiles\n',ifile,NbFiles);
    initialSound = DATA.rawSpeech{1,ifile};

    y = initialSound;
    y=y(:);%No noise
%     SigLength = length(y); %Length of the target signal
%     stdSig = std(y); % variance
%     
%     %%Frame by frame processing
% %     MFCC = zeros(nframe, 13);
%     MFCC = [];
%     n = 1;%Begining of a frame
%     m = frameLength;%End of a frame
%     iframe=1;
    [cepstra,aspectrum,pspectrum] = melfcc(y,Fs,D,...
        'lifterexp',0,...
        'nbands', 26,...
        'preemph',0,...
        'sumpower',1,...
        'maxfreq',8000);
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
