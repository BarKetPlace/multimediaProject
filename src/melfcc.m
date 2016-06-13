function [cepstra,aspectrum,pspectrum] = melfcc(y,En, sr,D, varargin)
%[cepstra,aspectrum,pspectrum] = melfcc(samples, sr[, opts ...])
%  Calculate Mel-frequency cepstral coefficients by:
%   - take the absolute value of the STFT
%   - warp to a Mel frequency scale
%   - take the DCT of the log-Mel-spectrum
%   - return the first <ncep> components
%  This version allows a lot of options to be controlled, as optional 
%  'name', value pairs from the 3rd argument on: (defaults in parens)
%    'wintime' (0.025): window length in sec
%    'hoptime' (0.010): step between successive windows in sec
%    'numcep'     (13): number of cepstra to return
%    'lifterexp' (0.6): exponent for liftering; 0 = none; < 0 = HTK sin lifter
%    'sumpower'    (1): 1 = sum abs(fft)^2; 0 = sum abs(fft)
%    'preemph'  (0.97): apply pre-emphasis filter [1 -preemph] (0 = none)
%    'dither'      (0): 1 = add offset to spectrum as if dither noise
%    'minfreq'     (0): lowest band edge of mel filters (Hz)
%    'maxfreq'  (4000): highest band edge of mel filters (Hz)
%    'nbands'     (40): number of warped spectral bands to use
%    'bwidth'    (1.0): width of aud spec filters relative to default
%    'dcttype'     (2): type of DCT used - 1 or 2 (or 3 for HTK or 4 for feac)
%    'fbtype'  ('mel'): frequency warp: 'mel','bark','htkmel','fcmel'
%    'usecmp'      (0): apply equal-loudness weighting and cube-root compr.
%    'modelorder'  (0): if > 0, fit a PLP model of this order
%    'broaden'     (0): flag to retain the (useless?) first and last bands
%    'useenergy'   (0): overwrite C0 with true log energy

% The following non-default values nearly duplicate Malcolm Slaney's mfcc
% (i.e. melfcc(d,16000,opts...) =~= log(10)*2*mfcc(d*(2^17),16000) )
%       'wintime': 0.016
%     'lifterexp': 0
%       'minfreq': 133.33
%       'maxfreq': 6855.6
%      'sumpower': 0
% The following non-default values nearly duplicate HTK's MFCC
% (i.e. melfcc(d,16000,opts...) =~= 2*htkmelfcc(:,[13,[1:12]])'
%  where HTK config has PREEMCOEF = 0.97, NUMCHANS = 20, CEPLIFTER = 22, 
%  NUMCEPS = 12, WINDOWSIZE = 250000.0, USEHAMMING = T, TARGETKIND = MFCC_0)
%     'lifterexp': -22
%        'nbands': 20
%       'maxfreq': 8000
%      'sumpower': 0
%        'fbtype': 'htkmel'
%       'dcttype': 3
% For more detail on reproducing other programs' outputs, see
% http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/mfccs.html
%
% 2005-04-19 dpwe@ee.columbia.edu after rastaplp.m.  
% Uses Mark Paskin's process_options.m from KPMtools
% $Header: /Users/dpwe/matlab/rastamat/RCS/melfcc.m,v 1.3 2012/09/03 14:01:26 dpwe Exp dpwe $


% 'nbands',26, 'useenergy',1,'maxfreq',8000,'hoptime'
% Parse out the optional arguments
[wintime, hoptime, numcep, lifterexp, sumpower, preemph, dither, ...
 minfreq, maxfreq, nbands, bwidth, dcttype, fbtype, usecmp, modelorder, ...
 broaden, useenergy] = ...
    process_options(varargin, 'wintime', 0.025, 'hoptime', 0.025, ...
          'numcep', 13, 'lifterexp', 0, 'sumpower', 1, 'preemph', 0, ...
	  'dither', 0, 'minfreq', 0, 'maxfreq', 8000, ...
	  'nbands', 26, 'bwidth', 1.0, 'dcttype', 2, ...
	  'fbtype', 'mel', 'usecmp', 0, 'modelorder', 0, ...
          'broaden', 0, 'useenergy', 0);

samples=y;

if preemph ~= 0
  samples = filter([1 -preemph], 1, samples);
end

% Compute FFT power spectrum
[pspectrum,logE] = powspec(samples, sr, wintime, hoptime, dither);
% pspectrum=power_spectrum(samples)
aspectrum = audspec(abs(pspectrum).^2, sr, nbands, fbtype, minfreq, maxfreq, sumpower, bwidth);

if (usecmp)
  % PLP-like weighting/compression
  aspectrum = postaud(aspectrum, maxfreq, fbtype, broaden);
end

%%%%%Denoise aspectrum (aspectrum is the mfccs before log and DCT
if ~isempty(D)
    [M,dsize]= size(D);
    
    Ey= aspectrum;
    
     mel_py = sum(abs(pspectrum).^2) ;
%     mel_px = sum(abs(pspectrumx).^2) ;
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
    
    %Silence frame    
    an = mel_py <= energythresh * ones(1, length(mel_py)) ;
    %Speech frames
    a = mel_py > energythresh * ones(1, length(mel_py)) ;

%     cd signalprocessing
%     [zhat]= getzhat(D, Ey, 60, En);
%     cd ..
    Ey_speech=Ey(:,a);
    En_speech= En(:,a);
    Ey_sil=Ey(:,an);    
    %Denoise silence
    Ey_sil= Ey_sil- mean(Ey(:,minI)); %*ones(1,nbframe);
    Ey_sil(Ey_sil<=1e-4)=1e-4;
    
    nbframe= size(Ey_speech,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        for iframe = 1:nbframe% Frame by frame processing
            fprintf('Frame %d/%d\n',iframe,nbframe);
            %
            zhat_frame= zeros(dsize,1);
            
            iproblem=1;
            ey=Ey_speech(:,iframe);
            en=En_speech(:,iframe);
            boundary= .02;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K= 15;
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
            
%             figure, plot(ey); hold on; plot(D*zhat_frame);
        end%   End Frame by frame processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aspectrum(:,a)=D*zhat;
     aspectrum(:,an)= Ey_sil;
%     figure, plot(Ey(:)); hold on; plot(aspectrum(:));
end


if modelorder > 0

  if (dcttype ~= 1) 
    disp(['warning: plp cepstra are implicitly dcttype 1 (not ', num2str(dcttype), ')']);
  end
  
  % LPC analysis 
  lpcas = dolpc(aspectrum, modelorder);

  % convert lpc to cepstra
  cepstra = lpc2cep(lpcas, numcep);

  % Return the auditory spectrum corresponding to the cepstra?
%  aspectrum = lpc2spec(lpcas, nbands);
  % else return the aspectrum that the cepstra are based on, prior to PLP

else
  
  % Convert to cepstra via DCT
  cepstra = spec2cep(aspectrum, numcep, dcttype);

end

cepstra = lifter(cepstra, lifterexp);

if useenergy
  cepstra(1,:) = logE;
end
  

