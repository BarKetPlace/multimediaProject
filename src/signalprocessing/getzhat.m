function [zhat] = getzhat(D,Ey,SNRtarget,En)
%[zhat] = getzhat(D,Ey,SNRtarget,En,mel_p)
%IN ::  D      -> Dictionnary (M x dsize matrix)
%       Ey      -> Observation (M x nbframe matrix)
%       SNRtarget -> snr target in dB (scalar)
%       En        -> Estimated noise on the feature (M x nbframe matrix)
%       mel_p     -> frame by frame power specrum (1 x nbframe vector)
%OUT::  zhat   -> estimate coefficients (dsize x nbframe matrix)
%The actual estimate of the feature Ey is D*zhat

[M,nbframe]= size(Ey);
dsize= size(D,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Frame by frame processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weight=1- mel_p./sum(mel_p);
% weight=ones(1,nbframe);


% En=(mean(En,2)*ones(1,nbframe)).*(ones(M,1)*weight);

% epsilon=sum((Ey-En).^2)/10^(SNRtarget/10);
if ~isempty(En)
    epsilon=max(sum(En.^2));
else
    epsilon=max(sum(Ey.^2)/10^(SNRtarget/10));
end
% extra= sum(En.^2);
% [Exhat, epsilon_tab, PrincipalCompNb, zhatstorage,SNR_Reconst]= ...
%                                     getEpsilon(nbframe,Ex, SNRtarget, D);
% Eyhat=zeros(size(Ey));
% zhat= zeros(dsize, nbframe);

for iframe = 1:nbframe
    iproblem=1;
    ey=Ey(:,iframe);
%     en=En(:,iframe);
    epsilon=max(sum(En.^2));
%     lambda=1e8;
    cvx_status='';
    while (~( strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) )
        cvx_begin quiet
            variables zhat_tmp(dsize) %pert(1)
            minimize( norm( zhat_tmp, 1 ))%+lambda*pert)% +100000*
            subject to
                D*zhat_tmp >= eps
                sum( (D*zhat_tmp-ey).^2) <= epsilon%+sum(en.^2)
        cvx_end
        

             fprintf('frame %d, Problem %d %s\n',iframe, iproblem,cvx_status);
        %     fprintf('Problem %d %s\n', iproblem,cvx_status);
        
        epsilon(1)=2^iproblem*epsilon(1);
        iproblem=iproblem+1;
    end
    %     zhatstorage(:,iframe)=zhat;
%      10*log10(sum(ex.^2)/sum((ex-D*zhat_tmp).^2))
    zhat(:,iframe)=zhat_tmp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Entire signal processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % En= mean(En,2)*ones(1,nbframe);
% 
% epsilon= sum((Ey-En).^2) /10^(SNRtarget/10);
% % weight=1- mel_p./sum(mel_p);
% % En=En.*(ones(M,1)*weight);
% 
% % % weight=ones(1,nbframe);
% % En=(En*ones(1,nbframe)).*(ones(M,1)*weight);
% iproblem=1;
%     cvx_status='';    
%     while (~( strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) )
%         cvx_begin
%         variables zhat(dsize,nbframe)
%         minimize( max( sum( abs(zhat) ) ) )
%         subject to
%         D*zhat >= eps
%         sum( (D*zhat-Ey+En).^2) <= epsilon%*ones(1,nbframe);
%         cvx_end
%         
%         epsilon=2^iproblem*epsilon;
%         iproblem=iproblem+1;
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('   \n');
% for iframe = 1:nbframe
%     fprintf('\b\b\b\b%02d%%\n',round(100*iframe/nbframe));
%     ex= Ex(:,iframe);
%
%     % snr_frame_mel_energy(iframe)= 10*log10(var(Ex)./var(En));
%     [ epsilon, frameSNR, zhat ] = getEpsilonFrame(ex, SNRtarget, D);
%
%     SNR_Reconst(iframe)=frameSNR;
%     epsilon_tab(iframe)=epsilon;
%     %Save results
%     [~,~,PrincipalCompNb(iframe)]= getPrincipalComp(zhat,.95);
%     %Save results
%     zhatstorage(:,iframe) = zhat;
%     Exhat(:,iframe) = D*zhat;
% end


end








% [M,dsize]=size(D); % Number of mfcc
% SparsityTarget = .5; %
% %Sparsity constrain
% [B] = lasso(D,ey);
% %Percentage of 0 coef in z
% perc=sum(B==0)/dsize;
% [~,idx]= min(abs(perc-SparsityTarget));
% zhat= B(:,idx);
% 
% Dtrunc=D(:,zhat~=0);
% 
% %Positivity constrain
% %We want to minimize ||ey-Dtrunc*zhat_lsq||^2 s.t. -Dtrunc*zhat_lsq <= -eps
% C=Dtrunc;
% d=ey;
% A=-Dtrunc;
% b=zeros(1,M)-eps;
% %No bound
% ub=inf(size(C,2),1);
% lb=-1*ub;
% %We choose the right algorithm that works for underdetermined system (less
% %equations than variable)
% options = optimoptions('lsqlin','Algorithm','active-set');
% 
% [zhat_lsq]= lsqlin(C,d,A,b,[],[],lb,ub,[],options);
% 
% % eyhat=Dtrunc*zhat_lsq;
% % tmp = find(eyhat<=0);
% % if any(tmp)
% %     warning('Non positive values');
% % end
% zhat(find(zhat~=0))=zhat_lsq;




