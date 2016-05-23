%% load signals and set parameters
close all
clc

    % get a signal from training set
load dataTrain.mat
Fs = 16000 ;                                                    % sampling frequency
x = DATA.rawSpeech{1, 325} ;                                    % get signal in DATA
x = x - mean(x) ;                                               % normalize
x = x / max(abs(x)) ;                                           %
[~, Ex, pspectrumx] = melfcc(x, Fs, [],'useenergy', 1);         % compute melfeatures

    % perform speech/silence recognition with threshold
energythresh = .2 ;                                             % threshold for speech/silence decision
mel_p = sum(pspectrumx) ;                                       
a = mel_p > energythresh * ones(1, length(mel_p)) ;
Ex= Ex(:,a);

    % load codebook of 6 bits
load Codebooks
D = Codebooks{1, 5} ;
[M, dsize] = size(D) ;
[D, D_par] = mapstd(D) ;

    % normalize Ex with training data parameters
mean_train = double(D_par.xmean) ;
std_train = double(D_par.xstd) ;
mean_train = repmat(mean_train, [1 length(Ex)]) ;
std_train = repmat(std_train, [1 length(Ex)]) ;
Ex = (Ex - mean_train) ./ std_train ; 

    % define SNR target;
SNRtarget = 20 ;            % dB

%% getEpsilon().m


    % set parameters
nframes = length(Ex) ;
epsilon = 0.1 ;
SNR = zeros(1, nframes) ;

% for i = 1:nframes
   
    ex = Ex(:, i) ;
    
    cvx_begin quiet
        variables zhat(dsize)
        minimize( norm( zhat, 1 ) )
        subject to
            D*zhat >= eps
            norm( D * zhat - ex, 2 ) <= epsilon
    cvx_end
    
    SNR(i) = 10*log10((norm(ex, 2).^2) / (norm(ex - D*zhat, 2).^2)) ;    
    
% end


