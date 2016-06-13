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
D=Codebooks{1,1};%
clear Codebooks;
[M, dsize]=size(D);


%Choose data
load ../dataTest.mat
% ISIGNAL= [round((length(DATA.rawSpeech)-1)*rand(1,20))+1];
% ISIGNAL= [10:20];
ISIGNAL=26;
sparsity=[];
t_=[];
En_=[];
En_model_=[];
ifig=1;

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
    %Extract mfcc
    [cepstrax,Ex,pspectrumx] = melfcc(x, [], Fs, []);
    [cepstray,Ey,pspectrumy] = melfcc(y, [], Fs, []);
    [cepstran,En_model,pspectrumn] = melfcc(n, [], Fs, []);
    cd signalprocessing
    fprintf('done.\n');
    
    %Compute power of each frame
    mel_py = sum(abs(pspectrumy).^2) ;
    mel_px = sum(abs(pspectrumx).^2) ;
    [~,minI]= min(mel_py);% Assumption:: should be zero

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      isolate filter bank energies that correspond to speech signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %find the energy threshold by
    i=2;
    while mel_py(i)<=2*mel_py(i-1)
        i=i+1;
    end
    energythresh=max(mel_py(1:i-1));      % threshold for speech/silence decision
    
    
%     energythresh=0; % Used only in plots
    
    
        
    %Silence frame
    an = mel_py <= energythresh * ones(1, length(mel_py)) ;
    % En_estimated=Ey(:,an);c
    En_silence= Ey(:,an);
    
    %Speech frames
    a = mel_py > energythresh * ones(1, length(mel_py)) ;
    
    speechmel_py=mel_py(a)/max(abs(mel_py(a)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % %Actual noise on the features
    En=Ey-Ex;
    
    %Store noises in case we process several signals
    En_=horzcat(En_,En);                     %Real noise 
    En_model_=horzcat(En_model_,En_model);   %Feature of the noise signal

%%
    nbframe=size(Ey(a),2);                      %Number of processed frame
    
    Ey_speech=Ey(:,a);
    En_speech=En(:,a);
    Ey_sil=Ey(:,an);
    En_sil=En(:,an);
    
    
    %Denoise silence
    Ey_sil= Ey_sil- mean(Ey(:,minI)); %*ones(1,nbframe);
    Ey_sil(Ey_sil<1e-4)=1e-4;
    %% Find the boundary epsilon
    %The epsilon boudary is easy to find:
    %We want epsilon such that ||Ex-Exhat||_2<= epsilon  and
    %||Ex||_2/||Ex-Exhat||_2 >= SNRtarget in dB
      
    
    SNRtarget=40;%dB
    %% Find zhat s.t. zhat= min_z ||ey- D*z||, subject to D*z>=0
    % zhat=getzhat(D,Ey,K,En);
    

    K_= [15];
    for iK=1:length(K_)
        K=K_(iK);
        boundary=.02;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        for iframe = 1:nbframe% Frame by frame processing
            fprintf('Frame %d/%d\n',iframe,nbframe);
            %
            zhat_frame= zeros(dsize,1);
            
            iproblem=1;
            ey=Ey_speech(:,iframe);
%             en=En_speech(:,iframe);
%             ex=Ex(:,iframe);
            
            %     [estimatedI, zhat_tmp, r] = getSupport(D,5,ey,en);
            I= [];                               %Set of indices step after step
            r= [];                               %set of residuals step after step
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         K= 7;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Initialization
            k= 1; I=[]; r=[]; Ip=[]; Iu=[]; index=[]; 
            zhat_tmp_storage=[];
            %       Kth greatest values of D'*ey
            [~,index]= sort(abs(D'*ey),'descend');
            I(:,k)= index(1:K);
            
            %       Estimate
            %Convex problem
            cvx_begin quiet
            variable zhat_tmp(K,1)
            minimize( norm(ey - D(:,I(:,k))*zhat_tmp )  )
            subject to
                D(:,I(:,k))*zhat_tmp >= eps
            cvx_end
            
            r(:,k)= ey- D(:,I(:,k))*zhat_tmp; %(pinv(D(:,I(:,k)))*y);
            zhat_tmp_storage(:,k)= zhat_tmp;
            
            while(1)
                k=k+1;
                [~,index]= sort(abs(D'*r(:,k-1)),'descend');
                Ip= index(1:K);
                
                Iu= union(I(:,k-1),Ip);
                %Convex problem
                cvx_begin quiet
                variable zhat_tmp(length(Iu),1)
                minimize( norm(ey - D(:,Iu)*zhat_tmp) )
                subject to
                    D(:,Iu)*zhat_tmp >= eps
                cvx_end
                
                xhat(Iu,1)= zhat_tmp;%pinv(D(:,Iu))*y;
                xhat(setdiff(1:dsize,Iu))=0;
                %    figure, plot(y); hold on; plot(D*xhat)
                
                [~,index]= sort(abs(xhat),'descend');
                I(:,k)= index(1:K);
                
                %Convex problem
                cvx_begin quiet
                variable zhat_tmp(K,1)
                minimize(  norm(ey - D(:,I(:,k))*zhat_tmp ) )
                subject to
                    D(:,I(:,k))*zhat_tmp >= eps
                cvx_end
                
                r(:,k)= ey - D(:,I(:,k))*zhat_tmp;%(pinv(D(:,I(:,k)))*ey);
                %    figure, plot(ey); hold on; plot(D(:,I(:,k))*zhat_tmp)
                zhat_tmp_storage(:,k)= zhat_tmp;
                
                
                if  ( norm(r(:,k)) <= boundary )
                    %         k= k - 1;
                    Ihat= I(:,k);
                    break;
                elseif ( norm(r(:,k),2) >= norm(r(:,k-1),2) ) || k>2*K
                    Ihat= I(:,k-1);
                    zhat_tmp= zhat_tmp_storage(:,k-1);
                    break;
                end
            end
            %         figure, plot(ey); hold on; plot(D(:,Ihat)*zhat_tmp); hold on; plot(ex);
            %         legend('ey', 'exhat','ex');
            zhat_frame(Ihat)= zhat_tmp;
            % Save result for the given frame
            zhat(:,iframe)=zhat_frame;
            
            
        end%   End Frame by frame processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %
        Exhat(:,a)= D*zhat;
        Exhat(:,an)= Ey_sil;
        
        
        
        %     PrincipalCompNb= zeros(1,nbframe);
        %     for iframe = 1:nbframe
        %         [~, ~, PrincipalCompNb(1,iframe) ] = getPrincipalComp(zhat(:,iframe), .9999);
        %     end
        %     sparsity=horzcat(sparsity,PrincipalCompNb);
        
%         err_ratio(iK)= norm(Ex(:) - Exhat(:),2)/norm(En(:),2);
        MSE(iK)= norm(Ex(:)-Exhat(:));
        
        
        
    end%  End isignal=ISIGNAL
end%  End variation of K
toc
%% PLOTS
figure(ifig),  ifig=ifig+1;
plot(K_,MSE); 
xlabel('Number of components');
ylabel('MSE');
title('MSE= ||Ex-Exhat||');
saveas(gcf,['K_vs_MSE_dsize' num2str(dsize) '.png']);

%%
plotMFCC(ifig,Ex,Ey,Exhat); ifig=ifig+1;
[snr_denoise, snr_mel_energy ]= plotSNR(ifig,Ex,Exhat,En); ifig=ifig+1;

figure(ifig), clf;  ifig=ifig+1;
histogram(sparsity)%,round(nbframe/2));
title({['Number of components representing .95% of energy in ' num2str(length(sparsity)) ' zhat vectors']});%['epsilon= ' num2str(epsilon)]});
% 
figure(ifig), clf;  ifig=ifig+1;
plot(mel_py,'LineWidth',2); hold on;
plot([1 length(mel_py)],energythresh*[1 1])%soundsc(x,Fs)

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
