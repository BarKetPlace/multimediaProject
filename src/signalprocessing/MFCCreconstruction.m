clear all
% close all
clc
% cd /home/antoine/Documents/multimediaProject/src/signalprocessing


isignal= 74;
SNR=5;
noise_path = '../../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose data
 load ../dataTrain.mat
x = DATA.rawSpeech{1,isignal};

clear DATA
% Choose codebook
load ../Codebooks
D=Codebooks{1,2};%
clear Codebooks;


[M, dsize]=size(D);
D=D./(ones(M,1)*sqrt(sum(D.^2)));

[Noise, Fs] = audioread([noise_path noise_file]);

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
%Noisy signal
y= x;

clear Noise
fprintf('MFCC extraction...');
cd ..
%Extract mfcc
[cepstrax,Ex,pspectrumx] = melfcc(x, Fs, [],'useenergy',1);
[cepstray,Ey,pspectrumy] = melfcc(y, Fs, [],'useenergy',1);
[cepstran,En_model,pspectrumn] = melfcc(n, Fs, [],'useenergy',1);
cd signalprocessing
fprintf('done.\n');
%extract Speech frame
mel_p = sum(pspectrumx) ;

% isolate filter bank energies that correspond to speech signal
energythresh = 20 ;                             % threshold for speech/silence decision

a = mel_p > energythresh * ones(1, length(mel_p)) ;
Ex= Ex(:,a);
Ey= Ey(:,a);
En_model= En_model(:,a);
% 
nbframe=size(Ex,2);


%Actual noise on the features
En= Ey- Ex;
%% Find the boundary epsilon
framelen=.025;
SNRtarget=25;%dB

[Exhat, epsilon_tab, PrincipalCompNb, zhatstorage]= ...
                                    getEpsilon(nbframes,Ex, SNRtarget, D)

err_ratio=  norm(Ex(:) - Exhat(:),2)^2/norm(En(:),2)^2;
%     err_ratio_1(ilambda)=var(Ex(:)-Exhat_1(:))/var(En(:));
    
%     ilambda=ilambda+1;


%Compute snr
snr_mel_energy= 10*log10(sum(Ex.^2)./sum(En.^2));
snr_mel_energy_model=10*log10(sum(Ex.^2)./sum(En_model.^2));
snr_denoise = 10*log10(sum(Ex.^2)./sum((Ex-Exhat).^2));

logPy=exp(cepstray(1,:));
% logPx=logPx+abs(min(logPx));
% fprintf('Overall sparsity= %02d%%\n',round(100*sum(zhatstorage(:)==0)/length(zhatstorage(:))));
%% PLOTS
ifig=1;

plotMFCC(ifig,Ex,Ey,Exhat);ifig=ifig+1;
plotSNR(ifig,Ex,Exhat,En);ifig=ifig+1;

figure, histogram(PrincipalCompNb);

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
