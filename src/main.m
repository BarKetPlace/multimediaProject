clear all
close all
clc

path_db = '../db/';
filename = 'female44.wav';

[initialSound,Fs] = audioread([path_db,filename]);
%reduce signal size (to reduce the computation time)
% x = initialSound(1:round(end/2));
x = initialSound;

N = length(x); %Length of the target signal
varx = var(x); % variance

frameLength_time = 20; %Frame length in ms
frameLength = 30/1000*Fs;%Frame length in samples
nframe = ceil(N/frameLength);


DFTlength = frameLength;
%Additive Noise
% y = x; %Here we do not add noise
varn = varx/1000000000;
%the noise follows N(0,varn)
noise = sqrt(varn)*randn(size(x));

%The observation
y = initialSound + noise;

%% Mel filter bank (Tiphanie report section C.2.3 - Filter Bank)
%Order of the filtering (number of triangles)
M = 26;
%Overlap between triangles, percentage of step in mel domain (between 0 and 1)
Overlap = 0;

[FilterBank] = MelCepstrumFilterBank(M, Fs, Overlap, DFTlength);

%Plot
figure, plot(sum(FilterBank))

% hold on;
% for k=1:M
%     plot(FilterBank(k,:),'-r'); hold on;
% end

xlabel('Samples'); ylabel('Amplitude');
title(['Mel filter, overlap in Mel domain = ' num2str(Overlap*100) '%']);

%% Frame by frame processing

MFCC_tmp = zeros(nframe,M);
MFCC = zeros(nframe, 13);
n = 1;%Begining of a frame
m = frameLength;%End of a frame
iframe=1;
while (m ~= N)
    yf = y(n:m);
    % DFT
    YF = abs(fft(yf)).^2;
    
    %Power spectrum PS is approximatly PSx + PSn ?
    PS = YF(1:round(DFTlength/2)+1);
       
%Mel filtering, 
%   For all triangles
    for i=1:M
        %Compute the mean of the power spectrum weighted by triangle
        MFCC_tmp(iframe,i) = 1/DFTlength*FilterBank(i,:)*PS;
    end
    
    
    %TODO, denoise MFCC_tmp using EM algorithm
    %...
    %MFCC_hat = ...
    
    %DCT
    MFCC_tmp(iframe,:) = dct(log10(MFCC_tmp(iframe,:)));
    MFCC(iframe,1:13) = MFCC_tmp(iframe,2:14);
    
    n = n + frameLength;
    m = min(N, m+frameLength);
    iframe=iframe + 1;
end

%% Plot MFCC
figure,
for i=1:nframe
    plot((i-1)*13+1:i*13,MFCC(i,:)); hold on;
end
title('frame by frame MFCCs');