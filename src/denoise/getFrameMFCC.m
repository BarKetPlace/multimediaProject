function [Ey, MFCC] = getFrameMFCC(frame,FilterBank)
%function [Ey, MFCC] = getFrameMFCC(frame,FilterBank)
%IN :: frame        -> waveform frame to process
%      FilterBank   -> Mel cepstrum filter bank
%OUT:: Ey           -> MFCC before log and DCT (used to denoise)
%      MFCC         -> Common MFCCs

    % DFT
    [M, DFTlength] = size(FilterBank);
    YF = abs(fft(frame)).^2;
    MFCC_tmp = zeros(1,M);
    %Power spectrum PS is approximatly PSx + PSn ?
    PS = YF(1:round(DFTlength));
    PS=PS(:); %Make it a column   
%Mel filtering, 
%   For all triangles
    for i=1:M
        %Compute the mean of the power spectrum weighted by triangle
        MFCC_tmp(i) = 1/DFTlength*FilterBank(i,:)*PS;
    end
    
    
    %TODO, denoise MFCC_tmp
    %...
    %MFCC_hat = ...
    
    %DCT
    Ey = MFCC_tmp;
    MFCC = dct(log10(MFCC_tmp));
    MFCC = MFCC(2:14);
end

