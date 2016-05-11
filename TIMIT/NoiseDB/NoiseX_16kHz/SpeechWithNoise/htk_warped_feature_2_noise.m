% This programme is to try different feature parameters
% This programme has to run seperately once fro training data and once for
% testing data

clear
close all


% If train flag is 1 then, training. If train flag is 0 then, testing.
train_flag=0;


% Noise specification 
NoiseType='babble';  % This specifies the noise type (babble, pink, volvo (car), m109 (tank), factory1, hfchannel)
SNR=0;   % SNR in dB

str1='~/RESEARCH3/NoiseDB/NoiseX_16kHz/';
str1=cat(2,str1,NoiseType);
str1=cat(2,str1,'_16kHz.wav');
Noise=wavread(str1);


% Specification about Fs, Window Length and Window Shift
Fs=16000;                        % Sampling Frequency
WinLen=0.032;                   % 32ms
WinShift=0.010;                % 10 ms


if (train_flag==1)
    
    % For Training 
    fin=fopen('~/HTK_My/HTK_TIMIT1/HTK_TIMIT/work/TRAIN.wavlist');      
    fout=fopen('~/HTK_My/HTK_TIMIT1/HTK_TIMIT/work/TRAIN.mfclist');
    
    NoOfFiles= 4620;
    
    Noise= 0 * Noise; 
    
else
    
    % For Testing 
    fin=fopen('~/HTK_My/HTK_TIMIT1/HTK_TIMIT/work/TEST.wavlist');      
    fout=fopen('~/HTK_My/HTK_TIMIT1/HTK_TIMIT/work/TEST.mfclist');
    
    NoOfFiles= 1680;
    
end


% Window length and shift

WinLenInSamples=WinLen*Fs;      % Window length in samples
WinShiftInSamples=WinShift*Fs;  % Window shift in samples



for i=1:1   % Here we need to correct; train -> 4620, test -> 1680
    i;
    filein=fgetl(fin);
    fileout=fgetl(fout);
    
    if mod(i,100)==0
        i
    end
    
    
    finput=fopen(filein,'r');
    input=fread(finput,'int16');
    input = input(513:end);
    fclose(finput);

    noise_file(input,SNR,Noise,train_flag)

%     O= make_wfcc_features_2_noise(input,WinLenInSamples,WinShiftInSamples,Fs,SNR,Noise,train_flag);
%     writehtk_new(fileout,O,32e-3,9);
%     plot(input)
%     i
%     pause
%     close;
    
end

fclose(fin);
fclose(fout);
