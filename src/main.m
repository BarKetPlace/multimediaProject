clear all

% db = '/home/antoine/Documents/multimediaProject/TIMIT/';
% filename = 'TRAIN/DR1/FCJF0/SI1657.WAV';
% [y, Fs] = audioread([db filename]);
% if system('test -f Dfull.mat') 
% load ../dataTrain.mat;
% Dfull = getDictionnary(dataTrain);
% save('Dfull.mat','Dfull');
% clear dataTrain
% else
% load Dfull.mat
% end
load dataTestNoisy5dB.mat
isignal = 34;
Fs = 16000;
y = DATA.rawSpeech{1,isignal};
load dataTest.mat
x = DATA.rawSpeech{1,isignal};
clear DATA
load Codebooks; 
% soundsc(y,Fs);
% 10*log10(var(x)/(var(y)-var(x)))

%%
D=[];
[cepstrax,aspectrumx,pspectrumx] = melfcc(x, Fs, D,...
        'lifterexp',0,...
        'nbands', 26,...
        'preemph',0,...
        'maxfreq',8000);
[cepstray,aspectrumy,pspectrumy] = melfcc(y, Fs, D,...
        'lifterexp',0,...
        'nbands', 26,...
        'preemph',0,...
        'maxfreq',8000);
M=size(aspectrumx,1);
%%Normalization of the mfcc
% aspectrumx=aspectrumx-ones(M,1)*mean(aspectrumx);
% aspectrumx=aspectrumx./(ones(M,1)*max(abs(aspectrumx)));
% 
% aspectrumy=aspectrumy-ones(M,1)*mean(aspectrumy);
% aspectrumy=aspectrumy./(ones(M,1)*max(abs(aspectrumy)));

%% Denoising
D = Codebooks{1,5};
dsize=size(D,2);
nframes=size(aspectrumy,2);
aspectrumy_d = zeros(size(aspectrumy));
zhatstorage=zeros(dsize,nframes);
lambda=.017;
fprintf('    \n');
for iframe = 1:nframes
    fprintf('\b\b\b\b%02d%%\n',round(iframe/nframes*100));
    ey = aspectrumy(:,iframe); %ey = ex + en
    cvx_begin quiet
        variables  zhat(dsize)
        minimize( norm( D * zhat - ey, 2 ) + lambda*norm( zhat, 1 ) )
        subject to
            D * zhat >= eps
    cvx_end
    
    zhatstorage(:,iframe) = zhat;
    aspectrumy_d(:,iframe) = D*zhat;
end
%%

frameidx=100;

figure(1), clf;
subplot(121);
    plot(aspectrumx(:,frameidx));hold on;
    plot(aspectrumy(:,frameidx)); hold on;
    plot(aspectrumy_d(:,frameidx));
    legend('clean', 'noisy','denoised');
    title('MFCC before log and dct');
subplot(122);
    stem(zhatstorage(:,frameidx));
    title('zhat');
    
figure(2),  clf;
    plot(cepstrax(:,frameidx)); hold on;
    plot(cepstray(:,frameidx));
    title('MFCC');
    legend('clean', 'noisy');
    
figure(3), clf;
    plot(pspectrumx(:,frameidx)); hold on;
    plot(pspectrumy(:,frameidx)); 
    title('Power spectrum of the signals');
    legend('clean', 'noisy');
%% plot 
frameidx=100;
figure(2), clf
plot(aspectrumx(:,frameidx)); hold on
plot(aspectrumy(:,frameidx));
legend('clean','noisy');
    








% 
% %%%
% %%Tiphanie's code
% 
% M=26;
% framelen_sec=.025;%seconds
% frameLengthSamples = framelen_sec*Fs;
% 
% 
% [FB, startFreq, centreFreq, endingFreq] = funct_Filterbanks(M,Fs,frameLengthSamples);
% 
% 
% x=x';
% SigLength=length(x);
% n = 1;%Begining of a frame
% m = frameLengthSamples;%End of a frame
% iframe=1;
% nframes=ceil(length(x)/frameLengthSamples);
% frames=zeros(nframes,frameLengthSamples);
% while (m ~= SigLength)
%     xf = x(n:m);
%     %MFCC extraction
%     framesx(iframe,:) = xf;
%     
%     n = n + frameLengthSamples;
%     m = min(SigLength, m+frameLengthSamples);
%     iframe=iframe + 1;
% end
% HamWin = hamming(frameLengthSamples)';
% dct_matrix = dctmtx(M);
% [~,mfccoutx] = funct_GetMfcc(framesx, FB, HamWin, dct_matrix);
% 
% y=y';
% SigLength=length(y);
% n = 1;%Begining of a frame
% m = frameLengthSamples;%End of a frame
% iframe=1;
% nframes=ceil(length(y)/frameLengthSamples);
% frames=zeros(nframes,frameLengthSamples);
% while (m ~= SigLength)
%     yf = y(n:m);
%     %MFCC extraction
%     framesy(iframe,:) = yf;
%     
%     n = n + frameLengthSamples;
%     m = min(SigLength, m+frameLengthSamples);
%     iframe=iframe + 1;
% end
% HamWin = hamming(frameLengthSamples)';
% dct_matrix = dctmtx(M);
% [~,mfccouty] = funct_GetMfcc(framesy, FB, HamWin, dct_matrix);
% mfccoutx=mfccoutx./(ones(M,1)*max(abs(mfccoutx)));
% mfccouty=mfccouty./(ones(M,1)*max(abs(mfccouty)));
% mfccoutx=mfccoutx-ones(M,1)*mean(mfccoutx);
% mfccouty=mfccouty-ones(M,1)*mean(mfccouty);
% figure(2), clf
% plot(mfccoutx(:,frameidx)); hold on;
% plot(mfccouty(:,frameidx));
% legend('clean','noisy');
