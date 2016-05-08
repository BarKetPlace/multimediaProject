function [MFCC] = getFrameMFCC(frame,FilterBank)
%function [MFCC] = getFrameMFCC(frame,FilterBank)


    % DFT
    [M, DFTlength] = size(FilterBank);
    YF = abs(fft(frame));
    MFCC_tmp = zeros(1,M);
    %Power spectrum PS is approximatly PSx + PSn 
    PS = YF(1:round(DFTlength));
       
%Mel filtering, 
%   For all triangles
    for i=1:M
        %Compute the mean of the power spectrum weighted by triangle
        MFCC_tmp(i) = 1/DFTlength*FilterBank(i,:)*PS;
    end
    
    
    %TODO, denoise MFCC_tmp using EM algorithm
    %...
    %MFCC_hat = ...
    
    %DCT
    MFCC = dct(log(MFCC_tmp));
    MFCC = MFCC(1:13);
end

