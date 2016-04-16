clear all
close all
clc

path_db = '../db/';
filename = 'female44.wav';

[x,Fs] = audioread([path_db,filename]);
x = x(1:round(end/2));%Part of signal to keep
N = length(x); %Length of the target signal
varx = var(x);

frameLength_time = 30; %Frame length in ms
frameLength = 30/1000*Fs;%Frame length in samples


%Additive Noise
% y = x; %Here we do not add noise
varn = varx/1000;
%the noise follows N(0,varn)
noise = sqrt(varn)*randn(size(x));

%The observation
y = x + noise;

%% Mel filter bank
 %Order of the filtering (number of triangles)
M = 26;
%Last frequency in mel domain
LastMelFreq = 2595*log10(1+Fs/2/700);
%Step in mel domain
delta = LastMelFreq/(M+1);
% The three following matrix represents the Mel filter bank of order M
MelKeyPoints = zeros(3,M);%Mel domain
FreqKeyPoints= zeros(3,M);%Normal frequency domain
SamplesKeyPoints= zeros(3,M);%Samples domain
%First row  : begining of filter
startPoint = 1;
%Second row : center
centerPoint = 2;
%Third row  : End
endPoint = 3;

%Computation of the key points
%Initialisation
MelKeyPoints(:,1) = [0;delta;delta+delta];
for i=2:M
    MelKeyPoints(:,i) = [MelKeyPoints(endPoint,i-1);... %StartPoint
                        MelKeyPoints(endPoint,i-1)+delta;...%centerPoint
                        MelKeyPoints(endPoint,i-1)+2*delta];%endPoint
end

%Define triangles
    
%% Frame by frame processing
n = 1;%Begining of a frame
m = frameLength;%End of a frame

while (m ~= N)
    yf = y(n:m);
    % DFT
    YF = fft(yf);
    YF = abs(YF(1:floor(end/2)));
    
    % Mel filtering
    
    n = n + frameLength;
    m = min(N, m+frameLength);
end