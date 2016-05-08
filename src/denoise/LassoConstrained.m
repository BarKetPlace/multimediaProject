clear all

load ../dataTest.mat
load Dfull.mat
%%
%Overlap between triangles, percentage of step in mel domain (between 0 and 1)
Overlap = .5;
Fs=16000;
frameLength_time = 30; %Frame length in ms
frameLength = frameLength_time/1000*Fs;%Frame length in samples
DFTlength = frameLength;

[FilterBank] = MelCepstrumFilterBank(Fs, Overlap, DFTlength);
[y] = dataTest.rawSpeech{1,1};
SigLength = length(y); %Length of the target signal
%%

%%Frame by frame processing
%     MFCC = zeros(nframe, 13);
Ey = [];
n = 1;%Begining of a frame
m = frameLength;%End of a frame
iframe=1;
while (m ~= SigLength)
    yf = y(n:m);
    %MFCC extraction
    Ey(:,iframe) = getFrameMFCC(yf,FilterBank);
    
    n = n + frameLength;
    m = min(SigLength, m+frameLength);
    iframe=iframe + 1;
end


%%
ex = Ey(:,1);
%restrain dic
D = Dfull(:,1:1000);

%Notation of the paper
%Initialisation
[B,FitInfo]= lasso(D,ex);
beta(:,1) = B(:,end); lambda(1)=FitInfo.Lambda(end);
A(:,1) = beta(:,1)~=0;

