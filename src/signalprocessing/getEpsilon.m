function [epsilon_tab]= getEpsilon(Ex, SNRtarget, D)
%INPUT:         Ey              targeted feature (M x nbframe matrix)
%               SNRtarget       wished SNR for each frame
%               D               Codebook (M x dsize matrix)
%               En              Noise (M x nbframe matrix)
%                
%OUTPUT
%               zhat            estimated coefficient
%               epsilon
%               Ehat           E estimation, E=D*zhat (M x nbframe matrix)
[M,nbframe]= size(Ex);
dsize= size(D,2);
epsilon_tab=zeros(1,nbframe);

% epsilon=sum(Ex.^2)/10^(SNRtarget/10);

% [Exhat, epsilon_tab, PrincipalCompNb, zhatstorage,SNR_Reconst]= ...
%                                     getEpsilon(nbframe,Ex, SNRtarget, D);
% Eyhat=zeros(size(Ex));
% zhat= zeros(dsize, nbframe);
%   en=En;
for iframe = 1:nbframe
    iproblem=1;
    ex=Ex(:,iframe);
  
    cvx_status='';
    epsilon=sum(ex.^2)/10^(SNRtarget/10);
    
    while (~( strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) )
        cvx_begin quiet
            variables zhat_tmp(dsize)
            minimize( norm( zhat_tmp, 1 ) )
            subject to
                D*zhat_tmp >= eps
                sum( ( D * zhat_tmp - ex).^2) <= epsilon
        cvx_end
        
        %     cvx_begin quiet
        %         variables zhat(dsize,nbframe)
        %         minimize( max( sum( abs(zhat) ) ) )
        %         subject to
        %             D*zhat >= eps
        %             sum( (D*zhat-Ex).^2) <= epsilon %*ones(1,nbframe);
        %     cvx_end
            fprintf('frame %d, Problem %d %s\n',iframe, iproblem,cvx_status);
        %     fprintf('Problem %d %s\n', iproblem,cvx_status);
        
        epsilon=2^iproblem*epsilon;
        iproblem=iproblem+1;
    end
    %     zhatstorage(:,iframe)=zhat;
%      10*log10(sum(ex.^2)/sum((ex-D*zhat_tmp).^2))
    epsilon_tab(iframe)=epsilon;
%     zhat(:,iframe)=zhat_tmp;

end
%     Eyhat=D*zhat;


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

