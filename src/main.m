clear all
close all
clc

path_db = '../db/';
filename = 'female44.wav';

[x,Fs] = audioread([path_db,filename]);
%reduce signal size (to faster the computation)
x = x(1:round(end/2));

N = length(x); %Length of the target signal
varx = var(x); % variance

frameLength_time = 20; %Frame length in ms
frameLength = 30/1000*Fs;%Frame length in samples

DFTlength = frameLength;
%Additive Noise
% y = x; %Here we do not add noise
varn = varx/100000000;
%the noise follows N(0,varn)
noise = sqrt(varn)*randn(size(x));

%The observation
y = x + noise;

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
% %Make the filter symetric
% symMelFilter = zeros(1,DFTlength);
% symMelFilter(1:round(DFTlength/2+1)) = MelFilter;
% symMelFilter(round(DFTlength/2+1)+1:end) = MelFilter(end:-1:4);
% figure, plot(symMelFilter);
%% Frame by frame processing
n = 1;%Begining of a frame
m = frameLength;%End of a frame
iframe=1;
while (m ~= N)
    yf = y(n:m);
    % DFT
    YF = abs(fft(yf)).^2;
    
    %Power spectrum PS is approximatly PSx + PSn
    PS = YF(1:round(DFTlength/2)+1);
    
    %TODO, denoise PS using EM algorithm
    %...
    %PShat = ...
    
    PShat = PS;
    
    
    
    % Mel filtering, 
    %For all triangles
    for i=1:M
        %Compute the mean of the power spectrum weighted by triangle
        MFCC_(iframe,i) = 1/DFTlength*FilterBank(i,:)*PShat;
    end
    
    %DCT
    MFCC_(iframe,:) = dct(log10(MFCC_(iframe,:)));
    MFCC(iframe,1:13) = MFCC_(iframe,2:14);
    n = n + frameLength;
    m = min(N, m+frameLength);
    iframe=iframe + 1;
end

%% Plot MFCC
figure,
for i=1:M
    plot((i-1)*13+1:i*13,MFCC(i,:)); hold on;
end
title('frame by frame MFCCs');