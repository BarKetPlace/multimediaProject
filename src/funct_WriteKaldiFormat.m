function [] = funct_WriteKaldiFormat(utt,feats,outputFolder,filename,N)
% Writes the features form Matlab structure to Kaldi format.
% Generates .ark and .scp files in indicated outputfolder;
%
% Inputs :
%- utt : names of the speech signals from data.utt
%- feats : features structure from data.feature
%- outputFolder : Folder in which to put generated files
%- filenames : string of the name for .ark file
%- N : nb of utterance
% Outputs : NONE
%
% Create Headers for ark writing
HEADER_MAT_TRAIN = cell(N,5);
HEADER_MAT_TRAIN(1,:) = {utt(1), size(feats{1},1),size(feats{1},2),1,size(feats{1},1)};
FEATURE_MAT_TRAIN = [];
FEATURE_MAT_TRAIN = [FEATURE_MAT_TRAIN; feats{1}];
for i=2:N
HEADER_MAT_TRAIN(i,:) = {utt(i), size(feats{i},1) , size(feats{1},2),...
HEADER_MAT_TRAIN{i-1,5}+1, HEADER_MAT_TRAIN{i-1,5} + size(feats{i},1)};
FEATURE_MAT_TRAIN = [FEATURE_MAT_TRAIN; feats{i}];
end
disp('HEADER_MAT_TRAIN and FEATURE_MAT_TRAIN computed');
% Create Mfcc directory
command = ['mkdir ' outputFolder];
[~,~]=system(command);
% Write for Kaldi
statArk = arkwrite(filename, HEADER_MAT_TRAIN, FEATURE_MAT_TRAIN);
% BE CAREFUL NAME OF FILES : .1. = 1 split
disp(['Status ark-writing:' num2str(statArk)]);
statScp = ark2scp(filename);
disp(['Status scp-writing:' num2str(statScp)]);
end
