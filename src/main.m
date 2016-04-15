clear all
close all
clc

path_db = '../db/';
filename = 'female44.wav';

[x,Fs] = audioread([path_db,filename]);
x = x(1:round(end/2));%Part of signal to keep
N = length(x); %Length of the target signal
varx = var(x);

flen_time = 30; %Frame length in ms
flen = 30/1000*Fs;%Frame length in samples


%Additive Noise
% y = x; %Here we do not add noise
varn = varx/1000;
noise = sqrt(varn)*randn(size(x));
%the noise follows N(0,varn)

y = x + noise;%Y represents the observation



%Frame by frame processing
n = 1;%Begining of a frame
m = flen;%End of a frame

% Mel filtering
f = [1:round(N/2)]*Fs/N;
MelCoef = 2595*log10(1+f/700);
MelCoef_samples = round(N*MelCoef/Fs);
    %Define triangles
    
while (m ~= N)
    yf = y(n:m);
    % DFT
    YF = fft(yf);
    YF = abs(YF(1:floor(end/2)));
    
    % Mel filtering
    
    n = n + flen;
    m = min(N, m+flen);
end