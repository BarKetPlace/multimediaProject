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

%%

load ../dataTest.mat

%Overlap between triangles, percentage of step in mel domain (between 0 and 1)
Overlap = .5;
Fs=16000;
frameLength_time = 20; %Frame length in ms
frameLength = frameLength_time/1000*Fs;%Frame length in samples
DFTlength = frameLength;

[FilterBank] = MelCepstrumFilterBank(Fs, Overlap, DFTlength);
[y] = dataTest.rawSpeech{1,100};
SigLength = length(y); %Length of the target signal

%%Frame by frame processing
%     MFCC = zeros(nframe, 13);
Ey = [];
n = 1;%Begining of a frame
m = frameLength;%End of a frame
iframe=1;
while (m ~= SigLength)
    yf = y(n:m);
    %MFCC extraction
    Ey(:,iframe) = getFrameMFCC(yf,FilterBank);
    
    n = n + frameLength;
    m = min(SigLength, m+frameLength);
    iframe=iframe + 1;
end
%% 
close all
zhat=[];
epsilon = .001;
for iframe = 1:10
clear zhat
%Randomly choose a frame to try
ey = Ey(:,ceil((size(Ey,2))*rand));

%restrain dic
iDictio=5;
load Codebooks.mat;
D=Codebooks{1,iDictio};
dsize = size(D,2);
M = size(D,1);
lambda=.01;

cvx_begin
    variable zhat(dsize)
    minimize( norm( D * zhat - ey, 2 ) + lambda*norm( zhat, 1 )  )
    subject to
        D * zhat >= 0
cvx_end

% figure(iframe), stem(zhat)
figure(iframe), clf;
subplot(121)
plot(ey,'LineWidth',2); hold on; plot(D*zhat);
title(['Sparsity ' num2str(round(sum(zhat~=0)/dsize*100)) '%']);
legend('Target','Estimated');
subplot(122)
    stem(zhat);
end
%%
% [zhat] = getzhat(D,ey);



% %% Lasso (Sparsity constrain)
% [B, FitInfo] = lasso(D,ey);
% [~, I] = min(FitInfo.MSE);
% 
% 
% zhat_lasso=zeros(dsize,size(B,2));
% eyhat_lasso=zeros(26,size(B,2));
% MSE_lasso=zeros(1,size(B,2));
% l0_lasso=zeros(1,size(B,2));
% for ilambda = 1:size(B,2)
%     
%     zhat_lasso(:,ilambda) = B(:,ilambda);
%     eyhat_lasso(:,ilambda) = D*zhat_lasso(:,ilambda);
%     tmp = eyhat_lasso(:,ilambda)<=0;
%     if any(tmp)
%         eyhat_lasso(:,ilambda) = eyhat_lasso(:,ilambda) + abs(min(eyhat_lasso(:,ilambda))); 
%     end
%     MSE_lasso(ilambda) = sum((ey(:)-eyhat_lasso(:,ilambda)).^2);
%     l0_lasso(ilambda) = sum(zhat_lasso(:,ilambda)==0);
% end
% 
% %%Least square (Inequality)
% 
% C=D;d=ey;A=-D;b=zeros(1,26);
% [zhat_lsq_tmp,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
% 
% 
% tot = sum(zhat_lsq_tmp.^2);
% [~,tmp] = sort(zhat_lsq_tmp.^2);
% perc=[floor(zhat_lsq_tmp(tmp(end)).^2/tot*100)];
% zhat_lsq = zeros(dsize,length(zhat_lsq_tmp));
% 
% eyhat_lsq = zeros(26,length(zhat_lsq_tmp));
% MSE_lsq = zeros(1,length(zhat_lsq_tmp));
% l0_lsq = zeros(1,length(zhat_lsq_tmp));
% for index=1:length(zhat_lsq_tmp)%Try several energy
% 
%     %Find the perc% 
% %     index=0;
% %     while sum(zhat_lsq_tmp(tmp(end-index:end)).^2)<=perc(index)/100*tot
% %         index=index+1;
% %     end
% 
%     zhat_lsq(tmp(length(zhat_lsq_tmp)-index+1:end),index)=zhat_lsq_tmp(tmp(length(zhat_lsq_tmp)-index+1:end));
%     eyhat_lsq(:,index) = D*zhat_lsq(:,index);
%     MSE_lsq(index) = sum((ey(:)-eyhat_lsq(:,index)).^2);
%     l0_lsq(index) = sum(zhat_lsq(:,index)==0);
% end
% %%PLOTS
% % close all
% figure(1), clf
% plot(l0_lsq/dsize,MSE_lsq); hold on;
% plot(l0_lasso/dsize,MSE_lasso);
% 
% xlabel('Sparsity (percentage of coef = 0 in zhat)');
% ylabel('MSE');
% legend('Least square','Lasso');
% 
% nplot=50;
% figure(2), clf
% plot(ey,'LineWidth',2); hold on;
% % plot(D*zhat_lsq_tmp); hold on;
% plot(eyhat_lsq(:, nplot)); hold on;
% plot(eyhat_lasso(:, nplot));
% legend('Target',...
%     ['LSQ algorithm, MSE:: ' num2str(MSE_lsq(nplot)) ' Sparsity:: ' num2str(floor(l0_lsq(nplot)/dsize*100))],...
%     ['LASSO algorithm, MSE:: ' num2str(MSE_lasso(nplot)) ' Sparsity:: ' num2str(floor(l0_lasso(nplot)/dsize*100))]...
%     );
