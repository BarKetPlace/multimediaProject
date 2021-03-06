function [] = loadDB(listfile_path, opt)
%function [] = loadDB(listfile_path, opt)
%   Usage: loadDB("data_train.list"/"data_test.list"/"data_dev.list")
%     to create "data_train.list","data_test.list" and "data_dev.list"
%     call ./listfiles.sh (from the matlab folder)


if strcmp(listfile_path,'data_train.list')
    DATA_filename = 'dataTrain.mat';
elseif strcmp(listfile_path,'data_test.list')
    DATA_filename = 'dataTest.mat';
elseif strcmp(listfile_path,'data_dev.list')
    DATA_filename = 'dataDev.mat';
else 
    fprintf('Usage: loadDB("data_train.list"/"data_test.list"/"data_dev.list")\n');
    return;
end

if nargin<2
    opt.nbSignal = 5000;
end

%clean in case it already exists


fprintf(['Creating ' DATA_filename '...']);
%Initialization
DATA = struct('utt',[],'rawSpeech',[],'speechframes',[],'mfcc',[],'part1',[],'ourFeature',[]);
fidBatch = fopen(listfile_path);
tlineBatch = fgetl(fidBatch);
i = 0;
while ischar(tlineBatch) && i < opt.nbSignal
%
% counter
i = i + 1;
% % Phone 
% PHN_file=[tlineBatch(1:end-3) 'PHN']; 
% PHNcell= importdata(PHN_file,' ');
% for iphone= 1:length(PHNcell)
%     splitted= strsplit(PHNcell{iphone,1}, {' '});
%     
% end
% Name of utterance (for Kaldi tool)
C = strsplit(tlineBatch,{'/','.'}) ;

lenC = length(C);
% dataTrain.utt{i} = [C{9} ' ' C{10} ' ' C{11} ' ' C{12}] ;
% Name with train drX is blocking in Kaldi
% Trial without it:
DATA.utt{i} = [C{lenC-2} '_' C{lenC-1}] ;
% dataTrain.utt{i} = tlineTrain;
% Read Timit
[y,fs,wrd,phn,ffx] =readsph(tlineBatch,['s','t','w'],-1);
phnmat_rng= fs*cell2mat(phn(:,1));
DATA.speechframes{i}= round(phnmat_rng(2,1):phnmat_rng(size(phnmat_rng,1),1));
%Conditioning
    y=y-mean(y);
    y=y/max(abs(y));
% Add to structure
DATA.rawSpeech{i} = y' ;


% raw speech in row vector
% Update line
tlineBatch = fgetl(fidBatch);
end

save(DATA_filename,'DATA','-v7.3') ;
%     
%     system('./listfiles.sh');%Create the list of files, should be modified according to system architecture
%     pathes{1} = 'data_train.list';%List of training files
%     pathes{2} = 'data_test.list';%List of testing files
%     pathes{3} = 'data_dev.list';%List of dev files
%     fprintf('Preparation of TIMIT database...');
%     funct_SpeechTimit(pathes);%Three .mat files are saved in the current folder
%     fprintf('done.\n');
fprintf('done.\n');
end


