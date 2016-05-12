%% Build dictionnaries from feature vectors before log and DCT
close all
clc

    % Load training data
load ../dataTrain.mat

    % set parameters
Fs = 16000 ;                                     % Sampling frequency of data
NbFiles = length(DATA.utt) ;                     % Get number of wav files
iter = 1 ;                                       % iteration in Dict
energythresh = 0.2 ;                             % threshold for speech/silence decision

    % Set dictionnary size
Numcol = 382622 ;                       
M = 26 ;
Dict = zeros(M, Numcol) ;               % column size computed on already computed MFCC structures (MUST INCLUDE ALGORITHM
                                        % IF CHANGE ON DATA SET!!)                       
          
for ifile = 1:NbFiles       % For all wav files in data

        % Get speech signal
    y = DATA.rawSpeech{1,ifile} ;
    y=y(:) ;                            % Put it in one single column
    
        % get mel features
    cd ../
    [~, mel_e, mel_p] = melfcc(y, Fs, []) ;
    cd signalprocessing/
    
        % get energy per frame
    mel_p = sum(mel_p) ;
    
        % isolate filter bank energies that correspond to speech signal
    a = mel_p > energythresh * ones(1, length(mel_p)) ;
    fbe_speech = mel_e(:, a) ;
    
        % store them in dictionnary
    Dict(:, iter:iter+length(fbe_speech)-1) = fbe_speech ;
    
    iter = iter + length(fbe_speech) ;
        
end

save('../Dictionnary.mat', 'Dict') ;

%% Build codebooks using kmeans algorithm
close all
clc

    % load data
load ../Dictionnary.mat

    % initialize parameters
bits = 2:8 ;
Codebooks = cell(1, length(bits)) ;

    % compute codebooks from clustering of Dictionnary
for i = 1 : length(bits)        
    [~, Cb] = kmeans(Dict', 2^(bits(i))) ;          % call built-in kmeans function
    Codebooks(1, i) =  {Cb'} ;                        % store result in cell        
end
    
    % save result
save('../Codebooks.mat', 'Codebooks') ;

%% Compare project k-means function and built-in function
close all
clc

        %%% INPUT %%%
    % data matrix
row = 2 ;                           % number of variables
column = 20000 ;                    % number of points in Data set
% DATA = rand(row, column) ;          % normally distributed
DATA = randn(row, column) ;         % uniformly distributed

    % number of cluster
num_clust = 64 ;

    % run project function
tic
[idxs, C] = kmeans_project(DATA, num_clust) ;
toc

    % compare with matlab algorithm
tic
[idxs_matlab, C_matlab] = kmeans(DATA', num_clust) ;
toc
idxs_matlab = idxs_matlab' ;
C_matlab = C_matlab' ;

    % plot result for 2 dimension data
plot(DATA(1, :), DATA(2, :), 'o') ;
hold on
plot(C(1, :), C(2, :), 'x') ;
plot(C_matlab(1, :), C_matlab(2, :), '+') ;
legend('DATA', 'Own clustering', 'Matlab clustering') ;
