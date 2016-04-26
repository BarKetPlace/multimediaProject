function [ FilterBank ] = MelCepstrumFilterBank(Fs, Overlap, DFTlength)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M = 26;

%Last frequency in mel domain
LastMelFreq = 2595*log10(1+Fs/2/700);
%Step in mel domain
delta = LastMelFreq/(M+1);
%Overlap in Mel domain

% The three following matrix represents the Mel filter bank of order M
MelKeyPoints = zeros(3,M);%Mel domain
FreqKeyPoints= zeros(3,M);%Normal frequency domain
SamplesKeyPoints= zeros(3,M);%Samples domain
%First row  : begining of triangle
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
    MelKeyPoints(:,i) = [MelKeyPoints(endPoint,i-1)-Overlap*delta;...%StartPoint
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

end

