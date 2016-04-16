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

%Last frequency in mel domain
LastMelFreq = 2595*log10(1+Fs/2/700);
%Step in mel domain
delta = LastMelFreq/(M+1);
%Overlap in Mel domain
Overlap = 0;%0.1*delta;
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

%Computation of the key points in Mel domain
%   Initialisation
MelKeyPoints(:,1) = [0;delta/2;delta];
%   
for i=2:M
    MelKeyPoints(:,i) = [MelKeyPoints(endPoint,i-1)-Overlap;...         %StartPoint
                        MelKeyPoints(endPoint,i-1)+delta/2;...  %centerPoint
                        MelKeyPoints(endPoint,i-1)+delta];      %endPoint
end
%Computation of the key points in normal frequency domain
FreqKeyPoints = 700*(10.^(MelKeyPoints(:,:)./2595) - 1);
%Computation of the key points in DFT scale
SamplesKeyPoints = round(FreqKeyPoints(:,:) * DFTlength/Fs);

%Compute the coefficient values, stored in FilterBank in the DFT scale
%domain
FilterBank = zeros(M, round(DFTlength/2 + 1));

%For all filters
for i = 1:M
    %triangle centered in c, starting in a and ending in b
    a = SamplesKeyPoints(startPoint,i);
    b = SamplesKeyPoints(endPoint,i);
    c = SamplesKeyPoints(centerPoint,i);
    
    %Temporary window to make the computation
    tmpWindow = zeros(1,round(DFTlength/2 + 1));
    % For all coeff in the DFT scale
    for k=0:round(DFTlength/2)
        if (a<=k && k<=c)
            tmpWindow(k+1) = 2/(c*b-a*c-b*a+a^2)*(k-a);
        elseif (c<k && k<=b)
            tmpWindow(k+1) = 2/(a*c-a*b-c*b+b^2)*(c-k) +2/(b-a);
        end
    end
    FilterBank(i,:) = tmpWindow;
end
MelFilter = zeros(1,round(DFTlength/2 + 1));

figure,
for i = 1:M
    plot(FilterBank(i,:),'-b'); hold on;
    MelFilter = MelFilter + FilterBank(i,:);
end

figure, plot(MelFilter);
% %Make the filter symetric
% symMelFilter = zeros(1,DFTlength);
% symMelFilter(1:round(DFTlength/2+1)) = MelFilter;
% symMelFilter(round(DFTlength/2+1)+1:end) = MelFilter(end:-1:4);
% figure, plot(symMelFilter);
%% Frame by frame processing
n = 1;%Begining of a frame
m = frameLength;%End of a frame

while (m ~= N)
    yf = y(n:m);
    % DFT
    YF = abs(fft(yf));
    YF = abs(YF(1:round(DFTlength/2)+1));
    
    % Mel filtering
    FilteredY = MelFilter.*YF';
    figure(10), plot(FilteredY);
    
    n = n + frameLength;
    m = min(N, m+frameLength);
end