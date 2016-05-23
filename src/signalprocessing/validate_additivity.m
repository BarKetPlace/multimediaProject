clear all
% cd /home/antoine/Documents/multimediaProject/src/signalprocessing


isignal= 73;
SNR=5;
noise_path = '../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose data
% load dataTest.mat
x = DATA.rawSpeech{1,isignal};

% clear DATA
%Choose codebook
% load ../Codebooks
% D=Codebooks{1,2};%
% clear Codebook;


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
y= x+ n;


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
energythresh = .2 ;                             % threshold for speech/silence decision

a = mel_p > energythresh * ones(1, length(mel_p)) ;
Ex= Ex(:,a);
Ey= Ey(:,a);
En_model= En_model(:,a);

nbframe=size(Ex,2);
%Actual noise on the features
En= Ey- Ex;
%% Computation of denoising
framelen=.025;
ifig=1;

% for lambda=.007
%For each frame

ilambda=1;
% lambda_tab= [.001:.005:.09 .01:.1:1 1.5:.5:4 ];% .01:.01:.1 .1:.1:1];
% lambda_tab= [.001:.005:.1];
% lambda_tab=[2:2:10];
lambda_tab= [.05];%:.02:.1];
% err_ratio=zeros(1,length(lambda_tab));
for lambda=lambda_tab
    display(['lambda= ' num2str(lambda)]);
    for iframe = 1:nbframe
        %     fprintf('%02d%%\n',round(iframe/nbframe*100));
        
        ey= Ey(:,iframe);%Set friendly variable names
        ex= Ex(:,iframe);
        en_model= En_model(:,iframe);
        en= En(:,iframe);
        
        %Frame SNR in time space
        timewindow= round(1+(iframe-1)*framelen*Fs:min(length(x),iframe*framelen*Fs));
        snr_time(iframe)= 10*log10( var(x(timewindow))/var(n(timewindow)) );
        
        % snr_frame_mel_energy(iframe)= 10*log10(var(Ex)./var(En));

        %First step denoising
        % With normalization
        cvx_begin quiet
            variables zhat(dsize)
            minimize( norm( zhat, 1 ) )
            subject to
                D*zhat >= eps
                norm( D * zhat - ey, 2 ) <= .1  %norm(En(:,iframe),2 )%.1
%                 zhat >= eps
%                 norm(zhat,2) <= 1
        cvx_end
        
%         nbprincipal=8;
%         [zhat,sortedzhat_idx]= getPrincipalComp(zhat, nbprincipal);
        
        %Save results
        zhatstorage(:,iframe) = zhat;
        Exhat(:,iframe) = D*zhat;
        
%         %Without normalization
%         cvx_begin quiet
%             variables zhat_1(dsize)
%             minimize( norm( D * zhat_1 - ey, 2 ) + lambda/10*norm( zhat_1, 1 ) )
%             subject to
%                 D * zhat_1 >= eps
% %                 norm(zhat,2) <= 1
%         cvx_end
%         zhatstorage_1(:,iframe) = zhat_1;
%         Exhat_1(:,iframe) = D*zhat_1;
        % %Keep the greatest part of the energy

        %
        %
        %
        % %Second part of denoising to ensure that D*zhat>=0
        % D_lowdim=D(:,sortedzhat_idx(1:nbprincipal));
        % cvx_begin quiet
        %     variables zhat_lowdim(nbprincipal)
        %     minimize( norm( D_lowdim * zhat_lowdim - ey, 2 ))% + lambda*norm( zhat, 1 ) )
        %     subject to
        %         D_lowdim * zhat_lowdim >= eps
        %         norm(zhat_lowdim,1)<=.8
        % %         norm(zhat_lowdim,1)>=1-eps
        % cvx_end
        % zhat(sortedzhat_idx(1:nbprincipal)) = zhat_lowdim;
    end
    %Compute estimation error
    err_ExExhat= sum( abs(Ex-Exhat).^2 );
    err_EyEx= sum( abs(Ey-Ex).^2 );
    err_additivity= abs(Ey -Ex -En_model)./abs(Ex +En_model);
    
    err_ratio(ilambda)=  norm(Ex(:) - Exhat(:),2)^2/norm(En(:),2)^2;
%     err_ratio_1(ilambda)=var(Ex(:)-Exhat_1(:))/var(En(:));
    
    ilambda=ilambda+1;
end

%Compute snr
snr_mel_energy= 10*log10(sum(Ex.^2)./sum(En.^2));
snr_mel_energy_model=10*log10(sum(Ex.^2)./sum(En_model.^2));
snr_denoise = 10*log10(sum(Ex.^2)./sum((Ex-Exhat).^2));

logPy=exp(cepstray(1,:));
% logPx=logPx+abs(min(logPx));
% fprintf('Overall sparsity= %02d%%\n',round(100*sum(zhatstorage(:)==0)/length(zhatstorage(:))));
%% PLOTS

%%%Plot denoising results
% 
ifig=1;
% framelist= [M*[1:nbframe 1:nbframe]; min(Ex(:))*ones(1,nbframe)  max(Ey(:))*ones(1,nbframe)];
% [~,idx]=sort(framelist(1,:));
% framelist=framelist(:,idx);


figure(ifig), clf; ifig=ifig+1;
stem(M*[1:nbframe], max(Ey(:))*ones(1,nbframe),'--k'); hold on;
plot(Ey(:),'LineWidth',2); hold on;
plot(Ex(:),'LineWidth',2); hold on;
plot(Exhat(:));hold on;
legend('Frames','Ey','Ex','Exhat');
xlabel('Mel space coefficent');
ylabel('Value f coefficient');
set(gca,'Xtick',[round(M/2):M:length(Ex(:))]);
set(gca,'XtickLabel',[1:nbframe]);
title('Denoising results');

% figure(ifig), clf; ifig=ifig+1;
% stem(M*[1:nbframe],max(Ey(:))*ones(1,nbframe),'--k'); hold on;
% plot(Ey(:));hold on;
% plot(Ex(:)+En_model(:));
% legend('Frames','Ey', 'Ex + En_{model}');
% title('Additivity check');

figure(ifig); clf; ifig=ifig+1;
% subplot(121)
% plot(err_ExExhat); hold on;
% plot(err_EyEx); hold on;
% legend('\Sigma|Exhat-Ex|^2','\Sigma|Ey-Ex|^2');
% title(['Lambda:: ' num2str(lambda)]);
% xlabel('Frame number'); 
% subplot(122);
plot(snr_denoise); hold on;
plot(snr_mel_energy); hold on;
legend('||Ex||_2/||Ex-Exhat||_2 in dB','||Ex||_2/||En||_2 in dB');
xlabel('Frame number'); ylabel('SNR in dB');
% [ptx,ptsy] = ginput(1);%(axis);%(choosepointfig);
% iframe = round(ptx);
title('SNR compare in Mel energy space');


figure(ifig); clf; ifig=ifig+1;
stem(dsize*[1:nbframe], max(zhatstorage(:))*ones(1,nbframe),'--k'); hold on;
stem(zhatstorage(:));
set(gca,'Xtick',[round(dsize/2):dsize:length(zhatstorage(:))]);
set(gca,'XtickLabel',[1:nbframe]);
title('zhat ');

if length(err_ratio)~=1
figure, plot(lambda_tab,err_ratio)
xlabel('lambda'); 
title({'Comparison btw denoise error and initial error';'ratio \sigma^2_{Ex-Exhat}/\sigma^2_{En}'});
end
