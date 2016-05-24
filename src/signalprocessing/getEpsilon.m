function [zhat, epsilon ]= getEpsilon(Ex, SNRtarget, D)

[M,nbframe]= size(Ex);                 
dsize= size(D,2);


epsilon= min( sum(Ex.^2)/10^(SNRtarget/10) );
% [Exhat, epsilon_tab, PrincipalCompNb, zhatstorage,SNR_Reconst]= ...
%                                     getEpsilon(nbframe,Ex, SNRtarget, D);
iproblem=1;
zhat= zeros(dsize, nbframe);
cvx_status='';
while (~( strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')) )
    cvx_begin quiet
        variables zhat(dsize,nbframe)
        minimize( max( sum( abs(zhat) ) ) )
        subject to
            D*zhat >= eps
            max( sum( (Ex - D*zhat).^2) )<= epsilon
    cvx_end
%     fprintf('Problem %s, epsilon= %f\n',cvx_status,epsilon);
    
    epsilon=2^iproblem*epsilon;
    iproblem=iproblem+1;
end
                            
                            

                            
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

