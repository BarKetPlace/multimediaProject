clear all
% close all
% clc
% cd /home/antoine/Documents/multimediaProject/src/signalprocessing



SNR=10;%dB             Wanted SNR of the noisy signal
noise_path = '../../TIMIT/NoiseDB/NoiseX_16kHz/';
noise_file = 'white_16kHz.wav';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Noise, Fs] = audioread([noise_path noise_file]);


% Choose codebook
load ../Codebooks.mat

EPSILON= [.001 .005 .01 .05 .1];
DSIZE= 1:3;
% for iepsilon=1:length(EPSILON)
    for iepsilon=3
    epsilon= EPSILON(iepsilon);
% for idsize=1:3
    for idsize=2
    
D= Codebooks{1,idsize};%
% clear Codebooks;
[M, dsize]=size(D);


%Choose data
load ../dataTest.mat
% ISIGNAL= [round((length(DATA.rawSpeech)-1)*rand(1,20))+1];
% ISIGNAL= [40:55];
ISIGNAL=47;
sparsity=[];
t_=[];
En_=[];
En_model_=[];
ifig=1;
eta=1e-4;

tic
for isignal=ISIGNAL
    %   get the signal from database
    x = DATA.rawSpeech{1,isignal};
    
    %   Keep the speech parts of the signal
%     x=x(DATA.speechframes{isignal});

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
    y= x+n;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     MFCC extraction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('MFCC extraction...');
    
    cd ..
    %Extract mfcc without denoising, 
    [~,Ex,pspectrumx] = melfcc(x, [], Fs, []);
    [~,Ey,pspectrumy] = melfcc(y, [], Fs, []);
    [~,En_model,pspectrumn] = melfcc(n, [], Fs, []);
    cd signalprocessing
    fprintf('done.\n');
    
%     Ex, Ey and E_model contains the Mel coefficient for the clean, noisy and noise signal.
%     We perform the denoising separately
    
    %Compute energy of each frame
    mel_py = sum(abs(pspectrumy).^2) ;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      isolate mel coefficients that correspond to speech signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %find the energy threshold
    i=2;
    while mel_py(i)<=2*mel_py(i-1)
        i=i+1;
    end
    energythresh=max(mel_py(1:i-1));      % threshold for speech/silence decision 
        
    %Silence frame
    an = mel_py <= energythresh * ones(1, length(mel_py)) ;

%     En_silence= Ey(:,an);
    
    %Speech frames
    a = mel_py > energythresh * ones(1, length(mel_py)) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % %Actual noise on the features
    En=Ey-Ex;
    
%     %Store noises in case we process several signals
%     En_=horzcat(En_,En);                     %Real noise 
%     En_model_=horzcat(En_model_,En_model);   %Feature of the noise signal

%%
    nbframe=size(Ey(:,a),2);                      %Number of speech frames
    
    
    Ey_speech=Ey(:,a);
    En_speech=En(:,a);
    
    Ey_sil= Ey(:,an);
    En_sil= En(:,an);
    Ex_sil= Ex(:,an);
    
    
    %Denoise silence by spectral substraction
    %Find the frame with the lowest energy
    [~,minI]= min(mel_py);
    
    Ey_sil= Ey_sil- mean(Ey(:,minI)); %*ones(1,nbframe);
    Ey_sil(Ey_sil<eta)=eta;% Normalize
    
    % Will receive the complete result
    Exhat= zeros(M, size(Ey,2));
    
    %Will only receive the processed speech frames
    Exhat_speech= zeros(M,nbframe);
    
    % Will receive the zhat values for all the processed speech frame
    zhat= zeros(dsize, nbframe);
        
    %% Find zhat s.t. zhat= min_z ||ey- D*z||, subject to D*z>=0
    % zhat=getzhat(D,Ey,K,En);
    

%     K_= [13];
%     for iK=1:length(K_)
%         K=K_(iK);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        for iframe = 1:nbframe% Frame by frame processing
%             fprintf('Frame %d/%d\n',iframe,nbframe);
            % Linear coefficient for the processed frame
            zhat_frame= zeros(dsize,1);
            
            % MFCCs of the processed frame
            ey=Ey_speech(:,iframe);
%%    VERSION 1 - SIMPLE CONVEX PROBLEM

            cvx_begin quiet
            variable zhat_frame(dsize,1)
            minimize( norm(zhat_frame,1) )
            subject to
                D*zhat_frame >= eta
                norm( ey-D*zhat_frame ) <= epsilon
            cvx_end
            
            %Unlikely but happen sometimes
            if (~( strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) )
                Exhat_speech(:,iframe)= Ey_speech(:,iframe);
                fprintf('frame %d:: Not solved\n',iframe);
            else% In case the problem is solved
                Exhat_speech(:,iframe)= D*zhat_frame;
                zhat(:,iframe)=zhat_frame;
            end
            
%%       VERSION 2 - IF number of Linear coefficient provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  K=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %Initialization
%             k= 1; I=[]; r=[]; Ip=[]; Iu=[]; index=[]; 
%             zhat_tmp_storage=[];
%             %       Kth greatest values of D'*ey
%             [~,index]= sort(abs(D'*ey),'descend');
%             I(:,k)= index(1:K);
%             
%             %       Estimate
%             %Convex problem
%             cvx_begin quiet
%             variable zhat_tmp(K,1)
%             minimize( norm(ey - D(:,I(:,k))*zhat_tmp )  )
%             subject to
%                 D(:,I(:,k))*zhat_tmp >= low_value
%             cvx_end
%             
%             r(:,k)= ey- D(:,I(:,k))*zhat_tmp; %(pinv(D(:,I(:,k)))*y);
%             zhat_tmp_storage(:,k)= zhat_tmp;
%             
%             while(1)
%                 k=k+1;
%                 [~,index]= sort(abs(D'*r(:,k-1)),'descend');
%                 Ip= index(1:K);
%                 
%                 Iu= union(I(:,k-1),Ip);
%                 %Convex problem
%                 cvx_begin quiet
%                 variable zhat_tmp(length(Iu),1)
%                 minimize( norm(ey - D(:,Iu)*zhat_tmp) )
%                 subject to
%                     D(:,Iu)*zhat_tmp >= low_value
%                 cvx_end
%                 zhat_frame(Iu,1)= zhat_tmp;%pinv(D(:,Iu))*y;
%                 
%                 [~,index]= sort(abs(zhat_frame),'descend');
%                 I(:,k)= index(1:K);
%                 
%                 %Convex problem
%                 cvx_begin quiet
%                 variable zhat_tmp(K,1)
%                 minimize(  norm(ey - D(:,I(:,k))*zhat_tmp ) )
%                 subject to
%                     D(:,I(:,k))*zhat_tmp >= low_value
%                 cvx_end
%                 
%                 r(:,k)= ey - D(:,I(:,k))*zhat_tmp;%(pinv(D(:,I(:,k)))*ey);
%                 %    figure, plot(ey); hold on; plot(D(:,I(:,k))*zhat_tmp)
%                 zhat_tmp_storage(:,k)= zhat_tmp;
%                 
%                 
%                 if  ( norm(r(:,k)) <= epsilon )
%                     %         k= k - 1;
%                     Ihat= I(:,k);
%                     break;
%                 elseif ( norm(r(:,k),2) >= norm(r(:,k-1),2) ) || k>2*K
%                     Ihat= I(:,k-1);
%                     zhat_tmp= zhat_tmp_storage(:,k-1);
%                     break;
%                 end
%             end
% 
%             zhat_frame(Ihat)= zhat_tmp;
%             % Save result for the given frame
%             zhat(:,iframe)=zhat_frame;
%             Exhat_speech= D*zhat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end%   End Frame by frame processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % 
        Exhat(:,a)= Exhat_speech;
        Exhat(:,an)= Ey_sil;
        
        
% Find the sparsity of the vector
        PrincipalCompNb= zeros(1,nbframe);
        for iframe = 1:nbframe
            [~, ~, PrincipalCompNb(1,iframe) ] = getPrincipalComp(zhat(:,iframe), .99);
        end
        sparsity=horzcat(sparsity,PrincipalCompNb);
        %Save sparsity results
        RES_mean(iepsilon,idsize)= mean(sparsity);
        RES_var(iepsilon,idsize)= var(sparsity);
        save('RES.mat','RES_mean','RES_var', 'EPSILON', 'DSIZE');
        
%         MSE= norm(Ex(:)-Exhat(:));
        
    end%  End isignal=ISIGNAL
 
toc
%% PLOTS
% figure(ifig),  ifig=ifig+1;
% plot(K_,MSE); 
% xlabel('Number of components');
% ylabel('MSE');
% title('MSE= ||Ex-Exhat||');
% saveas(gcf,['K_vs_MSE_dsize' num2str(dsize) '.png']);

%%
plotMFCC(ifig,Ex,Ey,Exhat); ifig=ifig+1;
[snr_denoise, snr_mel_energy ]= plotSNR(ifig,Ex,Exhat,En); ifig=ifig+1;

figure(ifig), clf;  ifig=ifig+1;
histogram(sparsity)%,round(nbframe/2));
title({['Nb of comp to get .99% of energy'];...
        ['dsize=' num2str(dsize) ', bound= ' num2str(epsilon)];...
        ['mean= ' num2str(mean(sparsity)) ', var= ' num2str(var(sparsity))]});
xlabel('Number of components');
saveas(gcf,['sparsity_dsize' num2str(dsize) '_bound' num2str(epsilon) '.png']);

end%End dsize
end
figure(ifig), clf;  ifig=ifig+1;
plot(mel_py,'LineWidth',2); hold on;
plot([1 length(mel_py)],energythresh*[1 1])%soundsc(x,Fs)
%% Kaldi results
WER= [74.5 75.3 76.2 77.9 79.0];
WER_BOUND= [.001 .01 .02 .05 .1];
NoProcessing= 80.5;
SilentProcessingPerf= 74.4;
CleanTestData= 56.3;
figure,
    plot(WER_BOUND, NoProcessing*ones(1,length(WER_BOUND)),'LineWidth',2,'Color','red'); hold on;
    plot(WER_BOUND, WER,'Marker','o','MarkerFaceColor','blue','Color','blue'); hold on;
    plot(WER_BOUND, SilentProcessingPerf*ones(1,length(WER_BOUND)),'LineWidth',2); hold on;
%     plot(WER_BOUND, CleanTestData*ones(1,length(WER_BOUND)),'LineWidth',2);
    
    legend('No processing','Silence + Speech frames','Only Silence frames');
    title({'Kaldi performances';'Noisy test data 10dB'});
    xlabel('Value of \epsilon');
    ylabel('WER (%)');

%%

figure, 
    plot(EPSILON,RES_mean(:,1),'LineWidth',2); hold on;
    plot(EPSILON,RES_mean(:,2),'LineWidth',2); hold on;
    plot(EPSILON,RES_mean(:,3),'LineWidth',2);
    legend('dsize= 32','dsize= 64','dsize= 128');
    ylabel('Mean of nb of comp in z');
    xlabel('epsilon value');
    title('Mean of nb of comp to get 99% of z');
    saveas(gcf,['Sparsity_mean.png']);
    
    figure, 
    plot(EPSILON,RES_var(:,1),'LineWidth',2); hold on;
    plot(EPSILON,RES_var(:,2),'LineWidth',2); hold on;
    plot(EPSILON,RES_var(:,3),'LineWidth',2);
    legend('dsize= 32','dsize= 64','dsize= 128');
    ylabel('Variance in number of components');
    xlabel('epsilon value');
    title('var of nb of comp to get 99% of z');
    saveas(gcf,['Sparsity_variance.png']);
% dsize= [32 64 128];
% epsilon= [.001 .01 .02];
% RES_mean= zeros(3,3); RES_var= zeros(3,3);
% RES_mean(1,:)= [14.3242 9.7896 8.1717]



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
