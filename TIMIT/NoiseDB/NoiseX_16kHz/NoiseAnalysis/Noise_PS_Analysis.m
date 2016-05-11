% This file checks how the noise spectrum varies

clear;
close all;

% Specification about Fs, Window Length and Window Shift
Fs=16000;                        % Sampling Frequency
WinLen=0.032;                   % 32ms
WinShift=0.010;                % 10 ms

% Window length and shift
WinLenInSamples=WinLen*Fs;      % Window length in samples
WinShiftInSamples=WinShift*Fs;  % Window shift in samples
Nshift=WinShiftInSamples;
HamWin=hamming(WinLenInSamples);




% Different Noise Types

% This specifies the noise type (babble, pink, volvo (car), m109 (tank), factory1, hfchannel)

noise_babble=wavread('~/RESEARCH3/NoiseDB/NoiseX_16kHz/babble_16kHz.wav');
noise_pink=wavread('~/RESEARCH3/NoiseDB/NoiseX_16kHz/pink_16kHz.wav');
noise_volvo=wavread('~/RESEARCH3/NoiseDB/NoiseX_16kHz/volvo_16kHz.wav');
noise_factory1=wavread('~/RESEARCH3/NoiseDB/NoiseX_16kHz/factory1_16kHz.wav');
noise_hfchannel=wavread('~/RESEARCH3/NoiseDB/NoiseX_16kHz/hfchannel_16kHz.wav');



% Finding the power spectrum
framecount=0;

start=1;
finish=WinLenInSamples;

omega=1:Fs/WinLenInSamples:Fs/2;

while (finish < (length(noise_babble)- WinLenInSamples))
    
    framecount=framecount+1
    
    % Windows of different noises
    sig_babble=noise_babble(start:finish);
    sig_pink=noise_pink(start:finish);
    sig_volvo=noise_volvo(start:finish);
    sig_factory1=noise_factory1(start:finish);
    sig_hfchannel=noise_hfchannel(start:finish);
    
    win_sig_babble = sig_babble .* HamWin; 
    win_sig_pink = sig_pink .* HamWin;
    win_sig_volvo = sig_volvo .* HamWin;
    win_sig_factory1 = sig_factory1 .* HamWin;
    win_sig_hfchannel = sig_hfchannel .* HamWin;
    
    PS_win_sig_babble = abs(fft(win_sig_babble)) .^2; 
    PS_win_sig_pink = abs(fft(win_sig_pink)) .^2; 
    PS_win_sig_volvo = abs(fft(win_sig_volvo)) .^2;
    PS_win_sig_factory1 = abs(fft(win_sig_factory1)) .^2;
    PS_win_sig_hfchannel = abs(fft(win_sig_hfchannel)) .^2;
    
%     PS_win_sig_babble = 10*log10( PS_win_sig_babble(1:(WinLenInSamples/2)) );
%     PS_win_sig_pink = 10*log10( PS_win_sig_pink(1:(WinLenInSamples/2)) );
%     PS_win_sig_volvo = 10*log10( PS_win_sig_volvo(1:(WinLenInSamples/2)) );
%     PS_win_sig_factory1 = 10*log10( PS_win_sig_factory1(1:(WinLenInSamples/2)) );
%     PS_win_sig_hfchannel = 10*log10( PS_win_sig_hfchannel(1:(WinLenInSamples/2)) );
    
    PS_win_sig_babble = PS_win_sig_babble(1:(WinLenInSamples/2)) ;
    PS_win_sig_pink =  PS_win_sig_pink(1:(WinLenInSamples/2)) ;
    PS_win_sig_volvo = PS_win_sig_volvo(1:(WinLenInSamples/2)) ;
    PS_win_sig_factory1 = PS_win_sig_factory1(1:(WinLenInSamples/2)) ;
    PS_win_sig_hfchannel = PS_win_sig_hfchannel(1:(WinLenInSamples/2)) ;
    
    
    
    % Plotting
    plot(omega,PS_win_sig_babble,omega,PS_win_sig_pink,omega,PS_win_sig_volvo,omega,PS_win_sig_factory1,...
        omega,PS_win_sig_hfchannel);
    legend('babble','pink','car','factory','hfchannel');
    
    pause;
    close;
    
       
    start=start+Nshift;
    finish=finish+Nshift;
        
end






