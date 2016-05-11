%% Build dictionnaries from feature vectors before log and DCT
clear all
    % Load training data
load ../dataTrain.mat
DATA = dataTrain ;
clear dataTrain;
    % Overlap between triangles, percentage of step in mel domain (between 0 and 1)
Overlap = .5 ;
Fs = 16000 ;

frameLength_time = 30 ;                  % Frame length in ms
frameLength = frameLength_time/1000*Fs ; % Frame length in samples
DFTlength = frameLength ;
[FilterBank] = MelCepstrumFilterBank(Fs, Overlap, DFTlength);   % Create filterbank
[M, ~] = size(FilterBank) ;
NbFiles = length(DATA.utt) ;                                    % Get number of wav files

    % Set dictionnary size
Numcol = 470747 ;
Dict = zeros(M, Numcol) ;               % column size computed on already computed MFCC structures (MUST INCLUDE ALGORITHM
                                        % IF CHANGE ON DATA SET!!)
k = 1 ;                                 % initialize number of column        

for ifile = 1:NbFiles       % For all wav files in data

        % Get speech signal and properties
    y = DATA.rawSpeech{1,ifile} ;
    y=y(:) ;                            % Put it in one single column
    SigLength = length(y);              % Get length
    varx = var(y);                      % Get Variance
    
        % Frame by frame processing of signal
    MFCC = [] ;
    n = 1 ;                      % Begining of a frame
    m = frameLength ;            % End of a frame
    iframe = 1 ;                 % Number of frame
    
    while (m ~= SigLength)       % Processes each frame
%         fprintf('Processing frame %d of signal %d/4620 | column %d/%d of dictionnary', iframe, ifile, k, Numcol) ;
        
            % Get frame of signal
        yf = y(n:m);
        
        [Ey, ~] = getFrameMFCC(yf,FilterBank);
        Dict(:,k) = Ey';
%         
%             % Compute Power spectrum
%         [M, DFTlength] = size(FilterBank);
%         YF = abs(fft(yf)).^2 ;               % absolute value of DFT
%         PS = YF(1:round(DFTlength)) ;
% 
%             %Mel filtering, 
%             %   For all triangles
%         for i=1:M
%             %Compute the mean of the power spectrum weighted by triangle
%             Dict(i, k) = 1/DFTlength*FilterBank(i,:)*PS;
%         end
%         
        
            % Update loop variables
        n = n + frameLength;
        m = min(SigLength, m+frameLength);
        iframe=iframe + 1;
        k = k + 1 ;
    end
end

%% Build k-means algorithm
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

    % number of iteration
Nb_iter = 100 ;

tic
        %%% FUNCTION %%%
    % 1) Initialization of centroids.
    % METHOD : choose randomly within the data (Forgy method)
C = DATA(:, randperm(column, num_clust)) ;          % get num_clust unique indices from all the possible indices

DATAmat = repmat(DATA, [1 1 num_clust]) ;           % initialize DATAmat matrix which will be used to compute distance
DATAmat = permute(DATAmat, [1 3 2]) ;               % permute dimension so that matrices fit

for j = 1 : Nb_iter
    % 2) For each point in the DATA, find indice of nearest centroid using euclidian
    % distance
Cmat = repmat(C, [1 1 column]) ;                    % add dimension to matrix to compute everything without for loop

dist = squeeze(sum((Cmat - DATAmat).^2)) ;          % compute euclidian distance between each point and each centroids

idxs = repmat(min(dist), [num_clust 1]) ;           % get the indice of the nearest centroid for each point
idxs = (1:num_clust) * (dist == idxs) ;             %

    % 3) Update the centroid of each cluster by computing the mean of all
    % the points contained in each cluster
for i = 1 : num_clust
   x = DATA(repmat(idxs == i, [row 1])) ;           % get all the points in cluster i
   x = reshape(x, [row length(x)/row]) ;            %
   
   C(:, i) = mean(x, 2) ;                           % computes the mean and stores it in centroid vector
end
end
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

%% Build codebooks using kmeans algorithm
close all
clc

    % load data
% load FeatDict.mat

    % initialize parameters
bits = 2:8 ;
Codebooks = cell(1, length(bits)) ;

    % compute codebooks from clustering of Dictionnary
for i = 1 : length(bits)        
    [~, Cb] = kmeans(Dict', 2^(bits(i))) ;          % call built-in kmeans function
    Codebooks(1, i) =  {Cb'} ;                        % store result in cell        
end







