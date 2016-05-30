clear all
% close all
% clc
% cd /home/antoine/Documents/multimediaProject/src/signalprocessing


isignal= 238;
SNR=5;
noise_path = '../../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Noise, Fs] = audioread([noise_path noise_file]);


% Choose codebook
load ../Codebooks
D=Codebooks{1,2};%
clear Codebooks;
[M, dsize]=size(D);
D=D./(ones(M,1)*sqrt(sum(D.^2)));


%Choose data
load ../dataTest.mat
ISIGNAL= [23:23];
ISIGNAL=63;
sparsity=[];
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
%% MFCC extraction
fprintf('MFCC extraction...');

cd ..
%Extract mfcc
[cepstrax,Ex,pspectrumx] = melfcc(x, Fs, [],'useenergy',1);
[cepstray,Ey,pspectrumy] = melfcc(y, Fs, [],'useenergy',1);
[cepstran,En_model,pspectrumn] = melfcc(n, Fs, [],'useenergy',1);
cd signalprocessing
fprintf('done.\n');

%Compute power of each frame
mel_p = sum(pspectrumy) ;
% figure, plot(mel_p); soundsc(x,Fs)
% % isolate filter bank energies that correspond to speech signal

%find the energy threshold by
i=2;
while mel_p(i)<=2*mel_p(i-1)
    i=i+1;
end
energythresh=max(mel_p(1:i-1));                             % threshold for speech/silence decision
% energythresh=0;
%Silence frame
an = mel_p <= energythresh * ones(1, length(mel_p)) ;
En_estimated=Ey(:,an);

%Speech frames
a = mel_p > energythresh * ones(1, length(mel_p)) ;
Ex= Ex(:,a);
Ey= Ey(:,a);
En_model= En_model(:,a);
speechmel_p=mel_p(a);
% %Actual noise on the features
En=Ey-Ex;

nbframe=size(Ex,2);

%% Find the boundary epsilon
%The epsilon boudary is easy to find: 
%We want epsilon such that ||Ex-Exhat||_2<= epsilon  ->
%||Ex||_2/||Ex-Exhat||_2 >= 20dB
SNRtarget=20;%dB
%It is computed in getzhat

%% Find zhat

zhat=getzhat(D,Ey,SNRtarget,En_estimated,speechmel_p);%En_estimated*ones(size(En))) ;
Exhat=D*zhat;


sigSNR= mean( 10*log10( sum(Ex.^2)./sum((Ex-Exhat).^2) ) );


PrincipalCompNb= zeros(1,nbframe);
for iframe = 1:nbframe
    [~, ~, PrincipalCompNb(1,iframe) ] = getPrincipalComp(zhat(:,iframe), .95);
end
sparsity=horzcat(sparsity,PrincipalCompNb);
% err_ratio= norm(Ex(:) - Exhat(:),2)^2/norm(En(:),2)^2;
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
% 
% %%
% % 
% % % figure(ifig), clf; ifig=ifig+1;
% % % plot(snr_mel_energy); hold on;
% % % plot(snr_time);
% % % 
% % % legend('mel energy domain','time domain');
% % % title('Frame by frame SNR');
% % % xlabel('Frame number');
% % % ylabel('SNRdB');
% % 
% % 
% % figure(ifig), clf; ifig=ifig+1;
% % plot(logPy); hold on;
% % plot(iframe*[1 1], [min(logPy) max(logPy)]);hold on;
% % % plot([1 nbframe], threshold*[1 1]);
% % xlabel('Frame number');
% % ylabel('sum(abs(X(f))^2)');
% % title('Power of each frame');
% % 
% % energy_snr=[snr_mel_energy;logPy];
% % 
% % [~,idx] = sort(snr_mel_energy);
% % % 
% % % 
% % % 
