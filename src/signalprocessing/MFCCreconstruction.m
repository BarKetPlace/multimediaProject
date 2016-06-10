clear all
% close all
% clc
% cd /home/antoine/Documents/multimediaProject/src/signalprocessing


isignal= 238;
SNR=10;
noise_path = '../../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Noise, Fs] = audioread([noise_path noise_file]);


% Choose codebook
load ../Codebooks
D=Codebooks{1,1};%
clear Codebooks;
[M, dsize]=size(D);
% D=D./(ones(M,1)*sqrt(sum(D.^2)));


%Choose data
load ../dataTest.mat
% ISIGNAL= [round((length(DATA.rawSpeech)-1)*rand(1,20))+1];
% ISIGNAL= [1:10];
ISIGNAL=46;
sparsity=[];
t_=[];
En_=[];
En_model_=[];
ifig=1;
tic
for isignal=ISIGNAL
x = DATA.rawSpeech{1,isignal};
n = Noise(1:length(x))';
%Conditioning of x
x=x-mean(x);
x=x/max(abs(x));
%extract the right noise length
noise_sig = Noise(1:length(x))';
%Uniformization of noise
UVnoise_sig = noise_sig/std(noise_sig);
UVnoise_sig = UVnoise_sig -mean(UVnoise_sig);
%future variance of noise depending on SNR
varn= (var(x)) / (10^(SNR/10));
%Amplification of noise
n = (varn^(.5))*UVnoise_sig;
% n=0;
% clear Noise

%Noisy signal
y= x+n;
% ffty=abs(fft(y)).^2;
% fftx=abs(fft(x)).^2;
% fftn=abs(fft(n)).^2;
% figure, plot(abs(fft(y)-fft(x)).^2, 'LineWidth',2); hold on; plot(fftn);
% figure, plot(fftx+fftn); hold on; plot(ffty);
%% MFCC extraction
fprintf('MFCC extraction...');

cd ..
%Extract mfcc
[cepstrax,Ex,pspectrumx] = melfcc(x, [], Fs, []);
[cepstray,Ey,pspectrumy] = melfcc(y, [], Fs, []);
[cepstran,En_model,pspectrumn] = melfcc(n, [], Fs, []);
cd signalprocessing
fprintf('done.\n');

%Compute power of each frame
mel_p = sum(abs(pspectrumy).^2) ;
% figure, plot(mel_p); soundsc(x,Fs)
% % isolate filter bank energies that correspond to speech signal

%find the energy threshold by
i=2;
while mel_p(i)<=2*mel_p(i-1)
    i=i+1;
end
energythresh=max(mel_p(1:i-1));                             % threshold for speech/silence decision
% energythresh=105;
%Silence frame
an = mel_p <= energythresh * ones(1, length(mel_p)) ;
% En_estimated=Ey(:,an);
En_silence= Ey(:,an);

% figure,histogram(En_silence(:))
%Speech frames
a = mel_p > energythresh * ones(1, length(mel_p)) ;
% Ex= Ex(:,a);
% Ey= Ey(:,a);
% En_model= En_model(:,a);
speechmel_p=mel_p(a)/max(abs(mel_p(a)));
% %Actual noise on the features
En=Ey-Ex;
En_=horzcat(En_,En);
En_model_=horzcat(En_model_,En_model);

nbframe=size(Ex,2); %Number of non-silence frame
%%
% t=mean(En_silence,2);

% cvx_begin quiet
%     variables t(M,1)%weight(1,nbframe)
%     minimize( norm(En - (t)*(1-speechmel_p),2) )
% cvx_end
% t_=horzcat(t_,t);

% t=[0.0282,0.0182,0.0619,0.0079,0.0110,0.1912,0.0865,0.0343, ...
%     0.0567,0.0178,0.0267,0.0184,0.0215,0.0195,0.0235,0.0286, ...
%     0.0324,0.0273,0.0272,0.0292,0.0276,0.0301,0.0283,0.0266, ...
%     0.0270,0.0247]';% Found by minimizing ||En-t*speechmel_p|| on the whole testing set

En_estimated= (median(En_silence,2))*(1-speechmel_p);
% En_estimated= En_model.*(ones(M,1)*speechmel_p);
% 10*log10( sum(En(:).^2)/sum( (En(:)-En_estimated(:)).^2) )
% figure(13), clf; plot(En(:)); hold on; plot(En_estimated(:));
%% Find the boundary epsilon
%The epsilon boudary is easy to find: 
%We want epsilon such that ||Ex-Exhat||_2<= epsilon  and
%||Ex||_2/||Ex-Exhat||_2 >= SNRtarget in dB
%It is computed in getzhat

SNRtarget=40;%dB
%% Find zhat

K_speech= 10;
K_sil= 10;
K= K_speech*ones(nbframe,1);
K(an)= K_sil;

zhat=getzhat(D,Ey,K,En);%En_estimated*ones(size(En))) ;
Exhat=D*zhat;


sigSNR= mean( 10*log10( sum(Ex.^2)./sum((Ex-Exhat).^2) ) );


PrincipalCompNb= zeros(1,nbframe);
for iframe = 1:nbframe
    [~, ~, PrincipalCompNb(1,iframe) ] = getPrincipalComp(zhat(:,iframe), .9999);
end
sparsity=horzcat(sparsity,PrincipalCompNb);
err_ratio(isignal)= norm(Ex(:) - Exhat(:),2)^2/norm(En(:),2)^2
%     err_ratio_1(ilambda)=var(Ex(:)-Exhat_1(:))/var(En(:));
    
%     ilambda=ilambda+1;


%Compute snr
% snr_mel_energy= 10*log10(sum(Ex.^2)./sum(En.^2));
% snr_mel_energy_model=10*log10(sum(Ex.^2)./sum(En_model.^2));
% snr_denoise = 10*log10(sum(Ex.^2)./sum((Ex-Exhat).^2));

% logPy=exp(cepstray(1,:));
% logPx=logPx+abs(min(logPx));
% fprintf('Overall sparsity= %02d%%\n',round(100*sum(zhatstorage(:)==0)/length(zhatstorage(:))));

end

toc
%% PLOTS

plotMFCC(ifig,Ex,Ey,Exhat); ifig=ifig+1;
[snr_denoise, snr_mel_energy ]= plotSNR(ifig,Ex,Exhat,En); ifig=ifig+1;
snr_mel_energy/snr_denoise

figure(ifig), clf;  ifig=ifig+1;
histogram(sparsity)%,round(nbframe/2));
title({['Number of components representing .95% of energy in ' num2str(length(sparsity)) ' zhat vectors']});%['epsilon= ' num2str(epsilon)]});
% 
figure(ifig), clf;  ifig=ifig+1;
plot(mel_p,'LineWidth',2); hold on;
plot([1 length(mel_p)],energythresh*[1 1])%soundsc(x,Fs)

% figure(ifig); clf; ifig=ifig+1;
% stem(dsize*[1:nbframe], max(zhatstorage(:))*ones(1,nbframe),'--k'); hold on;
% stem(zhatstorage(:));
% set(gca,'Xtick',[round(dsize/2):dsize:length(zhatstorage(:))]);
% set(gca,'XtickLabel',[1:nbframe]);
% title('zhat ');
% 
% % %% Frame by frame
% 
% iframe = round((nbframe-1)*rand())+1;
% iframe=9;
% 
% ey= Ey(:,iframe);%Set friendly variable
% ex= Ex(:,iframe);
% en_model= En_model(:,iframe);
% en= En(:,iframe);
% exhat=Exhat(:,iframe);
% if any(exhat<=0)
%     fprintf('Negative value\n');
% end
% 
% framefig=figure(ifig); clf; ifig=ifig+1;
% subplot(121);
% plot(ey,'LineWidth',2); hold on;
% plot(ex,'LineWidth',2); hold on
% plot(exhat);
% legend('Ey','Ex','Exhat');
% xlabel('Mel space coefficent');
% ylabel('Value f coefficient');
% subplot(122);
% stem(zhatstorage(:,iframe));
% 
% % 
% % %%
% % % 
% % % % figure(ifig), clf; ifig=ifig+1;
% % % % plot(snr_mel_energy); hold on;
% % % % plot(snr_time);
% % % % 
% % % % legend('mel energy domain','time domain');
% % % % title('Frame by frame SNR');
% % % % xlabel('Frame number');
% % % % ylabel('SNRdB');
% % % 
% % % 
% % % figure(ifig), clf; ifig=ifig+1;
% % % plot(logPy); hold on;
% % % plot(iframe*[1 1], [min(logPy) max(logPy)]);hold on;
% % % % plot([1 nbframe], threshold*[1 1]);
% % % xlabel('Frame number');
% % % ylabel('sum(abs(X(f))^2)');
% % % title('Power of each frame');
% % % 
% % % energy_snr=[snr_mel_energy;logPy];
% % % 
% % % [~,idx] = sort(snr_mel_energy);
% % % % 
% % % % 
% % % % 
