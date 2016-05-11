% This function is used to evaluate the features which is subsequently used
% by the HTK system

function O=make_wfcc_features_2_noise(input,WinLenInSamples,WinShiftInSamples,Fs,SNR,Noise,train_flag)



Nshift=WinShiftInSamples;

HamWin=hamming(WinLenInSamples)';



% MMFCC Feature parameters
M=26;                           % Number of triangular filters
Q=12;                           % Dimension of the feature vector (Default value)

M1=2;                           % Number of MFCC filters which need to be kept unchanged (i.e. K=700)
M2=M-M1;


% Warping constant
K1=700;                         % Default for MFCC
K2=900;                         % Warping constant for MMFCC


% Compression constant
a1=1;                           % Default for MFCC          
a2=0.1;                         % Compressing constant for MMFCC




s=input;


s=s-mean(s);

s= s / max(abs(s));


% ------------ Corrupting by noise ----------------


Noise_For_Input=Noise(1:length(s));

if (train_flag == 1)  % Then traning and no input noise (clean condition training)
    
    s_noisy=s;
    
else
    
    UV_Noise_For_Input=Noise_For_Input/std(Noise_For_Input);
    
    noise_var=(var(s) / (10^(SNR/10)));
    Noise_sample= (noise_var^(0.5)) * UV_Noise_For_Input;
    
    s_noisy=s+Noise_sample;
    
end


wavwrite(s,16000,'clean.wav');

wavwrite(s_noisy,16000,'clean_babble_snr10.wav');



