clear all

% db = '/home/antoine/Documents/multimediaProject/TIMIT/';
% filename = 'TRAIN/DR1/FCJF0/SI1657.WAV';
% [y, Fs] = audioread([db filename]);
load ../dataTrain.mat;
Dfull = getDictionnary(dataTrain);
clear dataTrain
%%

load ../dataTest.mat
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

[B,FitInfo] = lasso(D,ex);
[~, I] = min(FitInfo.MSE);
%%
close all

for delta = 4:10:size(B,2)
    figure(delta),
zhat = B(:,I+delta);
Eyhat = D*zhat;
%%Plot res

subplot(121); 
    plot(zhat);
   title({'Lambda'; num2str(FitInfo.Lambda(I+delta))});
subplot(122);
    plot(ex); hold on;
    plot(Eyhat);
title({'MSE'; num2str(FitInfo.MSE(I+delta))});
end