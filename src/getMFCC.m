function [MFCCcell] = getMFCC(dataTrain)
%% Mel filter bank (Tiphanie report section C.2.3 - Filter Bank)
%Order of the filtering (number of triangles)
M = 26;
%Overlap between triangles, percentage of step in mel domain (between 0 and 1)
Overlap = .5;
Fs = 16000;

frameLength_time = 30; %Frame length in ms
frameLength = frameLength_time/1000*Fs;%Frame length in samples
DFTlength = frameLength;
[FilterBank] = MelCepstrumFilterBank(M, Fs, Overlap, DFTlength);


NbFiles = length(dataTrain.utt);
fprintf('MFCC Extraction:     \n');
for ifile = 1:NbFiles
    fprintf('\b\b\b\b%02d%%\n',floor(ifile/NbFiles*100));
    initialSound = dataTrain.rawSpeech{1,ifile};
    %reduce signal size (to reduce the computation time)
    % x = initialSound(1:round(end/2));
    x = initialSound;
    y=x(:);%No noise
    SigLength = length(x); %Length of the target signal
    varx = var(x); % variance
    
    %%Frame by frame processing
%     MFCC = zeros(nframe, 13);
    MFCC = [];
    n = 1;%Begining of a frame
    m = frameLength;%End of a frame
    iframe=1;
    while (m ~= SigLength)
        yf = y(n:m);
        
        MFCC(1:13,iframe) = getFrameMFCC(yf,FilterBank);
        
        n = n + frameLength;
        m = min(SigLength, m+frameLength);
        iframe=iframe + 1;
    end
    
    dataTrain.mfcc{ifile} = MFCC;
end
%MFCC extracted, no need to carry the entire files

MFCCcell = dataTrain.mfcc;
end
