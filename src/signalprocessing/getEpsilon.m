function [Exhat, epsilon_tab, PrincipalCompNb, zhatstorage ]= getEpsilon(nbframes,Ex, SNRtarget, D)
for iframe = 1:nbframe
   
    ex= Ex(:,iframe);

    % snr_frame_mel_energy(iframe)= 10*log10(var(Ex)./var(En));
    [ epsilon, frameSNR, zhat ] = getEpsilonFrame(ex, SNRtarget, D);
    
    SNR(iframe)=frameSNR;
    epsilon_tab(iframe)=epsilon;
    %Save results
    [~,~,PrincipalCompNb(iframe)]= getPrincipalComp(zhat,.99);
    %Save results
    zhatstorage(:,iframe) = zhat;
    Exhat(:,iframe) = D*zhat;
end

end

