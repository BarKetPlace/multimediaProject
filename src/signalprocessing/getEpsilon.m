function [Exhat, epsilon_tab, PrincipalCompNb, zhatstorage, SNR_Reconst ]= ...
                                getEpsilon(nbframe,Ex, SNRtarget, D)

fprintf('   \n');
for iframe = 1:nbframe
    fprintf('\b\b\b\b%02d%%\n',round(100*iframe/nbframe));
    ex= Ex(:,iframe);

    % snr_frame_mel_energy(iframe)= 10*log10(var(Ex)./var(En));
    [ epsilon, frameSNR, zhat ] = getEpsilonFrame(ex, SNRtarget, D);
    
    SNR_Reconst(iframe)=frameSNR;
    epsilon_tab(iframe)=epsilon;
    %Save results
    [~,~,PrincipalCompNb(iframe)]= getPrincipalComp(zhat,.95);
    %Save results
    zhatstorage(:,iframe) = zhat;
    Exhat(:,iframe) = D*zhat;
end

end

