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

%% Build code books using dictionnary

    % yet to be coded