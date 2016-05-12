clear all
cd /home/antoine/Documents/multimediaProject/src/signalprocessing


isignal= 20;
SNR=5;
noise_path = '/home/antoine/Documents/multimediaProject/TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose data
load dataTest.mat
x = DATA.rawSpeech{1,isignal};
clear DATA
%Choose codebook
load ../Codebooks
D=Codebooks{1,2};%
clear Codebooks
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
y=x+n;




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
energythresh = 0.2 ;                             % threshold for speech/silence decision

a = mel_p > energythresh * ones(1, length(mel_p)) ;
Ex= Ex(:,a);
Ey= Ey(:,a);
En_model= En_model(:,a);

nbframe=size(Ex,2);
%Actual noise on the features
En= Ey-Ex;
%% Computation of SNR and denoising
framelen=.025;

%For each frame
for iframe = 1:nbframe
    fprintf('%02d%%\n',round(iframe/nbframe*100)); 

ey= Ey(:,iframe);%Set friendly variable
ex= Ex(:,iframe);
en_model= En_model(:,iframe);
en= En(:,iframe);




%SNR in Mel & time space
timewindow=round(1+(iframe-1)*framelen*Fs:min(length(x),iframe*framelen*Fs));
snr_time(iframe)=10*log10(var(x(timewindow))/var(n(timewindow)));

%Denoise
%denoising
[M dsize]=size(D);

lambda=.06;
    cvx_begin quiet
        variables  zhat(dsize)
        minimize( norm( D * zhat - ey, 2 ) + lambda*norm( zhat, 1 ) )
        subject to
            D * zhat >= eps
    cvx_end
zhatstorage(:,iframe) = zhat;
Exhat(:,iframe) = D*zhat;
% exhat=Exhat(:,iframe);

end
%
%Compute estimation error
err_ExExhat=sum( abs(Ex-Exhat).^2 );
err_EyExhat=sum( abs(Ey-Exhat).^2 );
err_additivity=mean(abs(Ey - Ex -En_model)./Ey);
%Compute snr
snr_mel_energy= 10*log10(var(Ex)./var(En));
snr_mel_energy_model=10*log10(var(Ex)./var(En_model));
snr_denoise = 10*log10(var(Ex)./var(Ex-Exhat));

logPy=exp(cepstray(1,:));
% logPx=logPx+abs(min(logPx));

%% PLOTS

%%%Plot denoising results
ifig=1;
figure(ifig), clf; ifig=ifig+1;
plot(err_additivity);

xlabel('Frame number');
% legend('Real noise', 'Model noise');
title('mean of relative error abs(ey -ex -en)/ey on a frame');

figure(ifig), clf; ifig=ifig+1;
plot(err_ExExhat,'LineWidth',2); hold on;
plot(err_EyExhat,'LineWidth',2); hold on;
% plot(exhat);
legend('Ex-Exhat','Ey-Exhat');
title(['Lambda:: ' num2str(lambda)]);
xlabel('Frame number'); 

figure(ifig), clf; ifig=ifig+1;
plot(snr_denoise); hold on;
plot(snr_mel_energy); hold on;
legend('var(Ex)/var(Ex-Exhat) in dB','var(Ex)/var(En) in dB');
xlabel('Frame number'); ylabel('SNR in dB');
title('SNR compare in Mel energy space');
%%
% 
% iframe = round((nbframe-1)*rand())+1;
% 
% ey= Ey(:,iframe);%Set friendly variable
% ex= Ex(:,iframe);
% en_model= En_model(:,iframe);
% en= En(:,iframe);
% exhat=Exhat(:,iframe);
% 
% figure(ifig), clf; ifig=ifig+1;
% plot(ey,'LineWidth',2); hold on;
% plot(ex,'LineWidth',2); hold on
% plot(exhat); hold on;
% % plot(en);
% legend('Ey','Ex','Exhat');
% xlabel('Mel space coefficent');
% ylabel('Value f coefficient');
% % title();
% 
% figure(ifig), clf, ifig=ifig+1;
% stem(zhatstorage(:,iframe));
% 
% % 
% % figure(ifig), clf; ifig=ifig+1;
% % plot(snr_mel_energy); hold on;
% % plot(snr_time);
% % 
% % legend('mel energy domain','time domain');
% % title('Frame by frame SNR');
% % xlabel('Frame number');
% % ylabel('SNRdB');
% 
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
% 
% 
% 
